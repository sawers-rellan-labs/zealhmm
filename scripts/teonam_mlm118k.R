#!/usr/bin/env Rscript
# =============================================================================
# TeoNAM STAM MLM (Q = 5 PCs fixed + K = Centered_IBS random, P3D) on the
# AUTHENTIC 118,838-SNP Chen 2019 GWAS panel (v5) — the faithful Fig 4C MLM.
#
# NO DENSIFICATION / INTERPOLATION. The 118K GWAS panel is only ~2.7% missing
# (Chen's MR < 0.2 SNPs), so the MLM runs DIRECTLY on it and PRESERVES the natural
# per-SNP scatter of Fig 4C. TASSEL handles per-site missingness (drops missing
# taxa per marker); K (Centered_IBS) uses observed allele frequencies. This is the
# whole point of the 118K panel — contrast scripts/teonam_mlm_interp.R, which
# step-interpolated the 73%-missing 51K MAP panel to a COMPLETE matrix and thus
# produced a block-terraced plateau Manhattan (a model demo, not the scatter).
#
# Q = 5 PCs from prcomp on a MEAN-IMPUTED copy of the dosage (PCA needs no NA);
# the MLM scan itself uses the raw ~2.7%-missing genotypes.
#
# Same TASSEL MLM invocation as teonam_mlm_interp.R (P3D, no compression).
# Stages: build | kinship | mlm | all   (default all)
# Inputs : data/teonam/teonam_gwas118k_dosage.rds  (markers x lines, NA=untyped)
#          data/teonam/markers_v5_gwas118k.tsv      (v5 chr/pos roster)
#          data/teonam/9250682/...phenotype...xlsx  (STAM)
# Outputs: data/teonam/tassel/geno_gwas_118k.hmp.txt, pheno_stam_pc_118k.txt,
#          kinship_cIBS_118k.txt, mlm_118k{1,2,3}.txt
# Run: Rscript scripts/teonam_mlm118k.R [all]
# =============================================================================
suppressMessages({
  library(data.table)
})
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)
source(file.path(ROOT, "scripts/logging.R"))
TASSEL <- "/Applications/TASSEL 5/run_pipeline.pl"
TDIR <- file.path(ROOT, "data/teonam/tassel")
OUT <- file.path(ROOT, "results/sim/teonam")
HMP <- file.path(TDIR, "geno_gwas_118k.hmp.txt")
PHENO <- file.path(TDIR, "pheno_stam_pc_118k.txt")
KIN <- file.path(TDIR, "kinship_cIBS_118k.txt")
MLMOUT <- file.path(TDIR, "mlm_118k")
stage <- if (length(commandArgs(TRUE)) >= 1) commandArgs(TRUE)[1] else "all"
tass <- function(a, mem = "-Xmx12g") {
  t <- system.time(rc <- system(sprintf('"%s" %s %s', TASSEL, mem, a)))
  list(sec = as.numeric(t["elapsed"]), rc = rc)
}
timings <- list()

# ---- build the 118K hapmap (NO interpolation) + phenotype with 5 PCs --------
if (stage %in% c("build", "all") || !file.exists(HMP)) {
  log_info("build 118K hapmap (raw, N at missing; no densification) ...")
  g118 <- readRDS("data/teonam/teonam_gwas118k_dosage.rds")
  dos <- g118$dos # markers x lines, NA = untyped
  mc <- fread("data/teonam/markers_v5_gwas118k.tsv")
  chr_by <- setNames(mc$chr_v5, mc$marker)
  pos_by <- setNames(mc$pos_v5, mc$marker)
  mk <- intersect(rownames(dos), mc$marker) # lifted markers only
  G <- dos[mk, , drop = FALSE] # markers x lines, integer 0/1/2, NA
  lines <- colnames(G)
  log_info(
    "  %d lifted markers x %d lines; missing %.2f%%",
    nrow(G), ncol(G), 100 * mean(is.na(G))
  )

  # recode 0/1/2 -> A/M/C, NA -> N (TASSEL handles N per site)
  cn <- matrix("N", nrow(G), ncol(G))
  cn[!is.na(G) & G == 0L] <- "A"
  cn[!is.na(G) & G == 1L] <- "M"
  cn[!is.na(G) & G == 2L] <- "C"
  hd <- data.table(
    `rs#` = mk, alleles = "A/C", chrom = as.integer(chr_by[mk]),
    pos = as.integer(pos_by[mk]), strand = "+", `assembly#` = NA, center = NA,
    protLSID = NA, assayLSID = NA, panelLSID = NA, QCcode = NA
  )
  hmp <- cbind(hd, as.data.table(cn))
  setnames(hmp, c(names(hd), lines))
  setorder(hmp, chrom, pos)
  fwrite(hmp, HMP, sep = "\t", quote = FALSE, na = "N")
  log_info("  wrote %s", basename(HMP))

  # Q = 5 PCs: prcomp needs complete data -> mean-impute a copy (PCA only).
  Xm <- t(G) # lines x markers
  mu <- colMeans(Xm, na.rm = TRUE)
  v <- colMeans(Xm^2, na.rm = TRUE) - mu^2
  keepc <- is.finite(v) & v > 0 # drop invariant/all-NA markers
  Xm <- Xm[, keepc, drop = FALSE]
  mu <- mu[keepc]
  na_ij <- which(is.na(Xm), arr.ind = TRUE)
  if (nrow(na_ij)) Xm[na_ij] <- mu[na_ij[, 2]] # mean-impute for PCA
  pc <- prcomp(Xm, center = TRUE, scale. = FALSE)$x[, 1:5]
  rm(Xm)
  gc()

  ph <- as.data.frame(readxl::read_excel("data/teonam/9250682/TeoNAM_1257RILs_22traits_phenotype_data.xlsx"))
  names(ph)[1] <- "line"
  stam <- suppressWarnings(as.numeric(ph$STAM))
  names(stam) <- ph$line
  keep <- intersect(rownames(pc), names(stam)[is.finite(stam)])
  pt <- data.table(Taxa = keep, STAM = stam[keep], as.data.table(pc[keep, ]))
  setnames(pt, c("Taxa", "STAM", paste0("PC", 1:5)))
  writeLines(c(
    "<Phenotype>", paste(c("taxa", "data", rep("covariate", 5)), collapse = "\t"),
    paste(c("Taxa", "STAM", paste0("PC", 1:5)), collapse = "\t")
  ), PHENO)
  fwrite(pt, PHENO, sep = "\t", quote = FALSE, col.names = FALSE, append = TRUE)
  log_info("  wrote %s (%d taxa)", basename(PHENO), nrow(pt))
}

if (stage %in% c("kinship", "all")) {
  log_info("kinship (Centered_IBS on the raw 118K panel) ...")
  r <- tass(sprintf('-importGuess "%s" -KinshipPlugin -method Centered_IBS -endPlugin -export "%s" -exportType SqrMatrix', HMP, KIN))
  log_info("  kinship rc=%d, %.1f s", r$rc, r$sec)
  timings[["kinship"]] <- r$sec
}
if (stage %in% c("mlm", "all")) {
  log_info("MLM (P3D, Q+K) on the raw 118K panel ...")
  r <- tass(sprintf('-fork1 -h "%s" -fork2 -r "%s" -fork3 -k "%s" -combine4 -input1 -input2 -intersect -combine5 -input4 -input3 -mlm -mlmVarCompEst P3D -mlmCompressionLevel None -export "%s" -runfork1 -runfork2 -runfork3', HMP, PHENO, KIN, MLMOUT))
  log_info("  MLM rc=%d, %.1f s; files: %s", r$rc, r$sec, paste(list.files(TDIR, pattern = "mlm_118k"), collapse = ", "))
  timings[["mlm"]] <- r$sec
}
if (length(timings)) {
  tt <- data.table(step = names(timings), sec = unlist(timings))
  fwrite(tt, file.path(OUT, "mlm_118k_timings.csv"))
  print(tt)
}
