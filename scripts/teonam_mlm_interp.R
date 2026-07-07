#!/usr/bin/env Rscript
# Correct GWAS model (MLM: Q=5 PCs fixed + K=Centered_IBS random, P3D) on the
# INTERPOLATED COMPLETE genotype matrix (teonam_map_v5_gwas, 51,004 markers x
# 1,237 lines, step-interpolated -> no missing). Complete matrix => P3D reuses one
# kinship eigendecomposition across all markers => tractable (the sparse map panel
# was not). Same 51,004 markers as the interpolated OLS scan (stam_gwas_scan_interpolated.csv),
# so OLS-vs-MLM lambda_GC is apples-to-apples. Duplicate-cM markers -> identical
# interpolated rows (terraces); kept so every marker plots at its real bp.
# ACCEPTED CAVEAT: interpolation => block-constant genotypes => plateau Manhattan
# (no scatter). This run is about the MODEL (OLS->MLM Q+K), not the scatter.
# Run: Rscript agent/teonam_mlm_interp.R <build|kinship|mlm|all>
suppressMessages({
  library(data.table)
  devtools::load_all("/Users/fvrodriguez/repos/nilhmm", quiet = TRUE)
})
source("R/simulate.R")
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)
TASSEL <- "/Applications/TASSEL 5/run_pipeline.pl"
TDIR <- file.path(ROOT, "data/teonam/tassel")
OUT <- file.path(ROOT, "results/sim/teonam")
HMP <- file.path(TDIR, "geno_gwas_interp.hmp.txt")
PHENO <- file.path(TDIR, "pheno_stam_pc_interp.txt")
KIN <- file.path(TDIR, "kinship_cIBS_interp.txt")
MLMOUT <- file.path(TDIR, "mlm_interp")
stage <- if (length(commandArgs(TRUE)) >= 1) commandArgs(TRUE)[1] else "all"
tass <- function(a, mem = "-Xmx10g") {
  t <- system.time(rc <- system(sprintf('"%s" %s %s', TASSEL, mem, a)))
  list(sec = as.numeric(t["elapsed"]), rc = rc)
}
timings <- list()

# ---- build the interpolated complete matrix -> hapmap ----------------------
if (stage %in% c("build", "all") || !file.exists(HMP)) {
  message("densify (step-interpolate) to complete union matrix ...")
  mc <- fread("data/teonam/markers_v5.tsv") # roster + v5 bp (map-neutral liftover)
  setnames(mc, "chr_v5", "chr")
  nat_cm <- fread("data/teonam/teonam_v5_native.tsv") # native est.map: cM for placed markers
  mc[, cm := nat_cm$cm[match(marker, nat_cm$marker)]] # native cM; NA where est.map didn't place it
  mc[, cm := {
    ok <- !is.na(cm)
    if (any(!ok) && sum(ok) >= 2L) {
      f <- .bp_to_cm_fun(data.table(bp = pos_v5[ok], cm = cm[ok])) # native Marey spline (Hyman, monotone)
      cm[!ok] <- f(pos_v5[!ok])
    }
    cm
  }, by = chr] # place est.map-unplaced markers on the NATIVE cM scale via its Marey spline
  cm_by <- setNames(mc$cm, mc$marker)
  pos_by <- setNames(mc$pos_v5, mc$marker)
  # FULL 51,004 GWAS target (teonam_map_v5_gwas). interpolate_genotype accepts a
  # duplicate-cM target (only obs$cm must be strictly increasing), so we target all
  # genotyped markers directly -- duplicate-cM twins come out as identical rows.
  gcols <- names(fread("data/teonam/TeoNAM_genotype_clean.csv", nrows = 0))[-(1:3)]
  GWAS_MK <- intersect(gcols, mc$marker)
  fams <- c(
    TIL01 = "W22TIL01_genotype.csv", TIL03 = "W22TIL03_genotype.csv",
    TIL11 = "W22TIL11_genotype.csv", TIL14 = "W22TIL14_genotype.csv", TIL25 = "W22TIL25_genotype.csv"
  )
  fam_data <- lapply(names(fams), function(fam) {
    g <- fread(file.path("data/teonam", fams[fam]))
    g <- g[!duplicated(g[[1]])]
    list(g = g, keys = paste0(fam, sub("^.*Line_", "", g[[1]])))
  })
  names(fam_data) <- names(fams)
  densify_chr <- function(ch) {
    mch <- mc[chr == ch]
    setorder(mch, cm)
    tgt <- mch[marker %in% GWAS_MK]
    setorder(tgt, cm)
    tgt_df <- data.frame(chr = ch, cm = tgt$cm)
    blocks <- lapply(names(fam_data), function(fam) {
      g <- fam_data[[fam]]$g
      mk <- intersect(names(g)[-(1:3)], mch$marker)
      obs <- data.frame(marker = mk, cm = cm_by[mk])
      obs <- obs[order(obs$cm), ]
      obs <- obs[!duplicated(obs$cm), ]
      geno <- t(as.matrix(g[, obs$marker, with = FALSE]))
      storage.mode(geno) <- "double"
      dn <- interpolate_genotype(geno, data.frame(chr = ch, cm = obs$cm), tgt_df, mode = "step")
      colnames(dn) <- fam_data[[fam]]$keys
      dn
    })
    D <- do.call(cbind, blocks)
    rownames(D) <- tgt$marker
    list(geno = D, markers = tgt[, .(marker, chr, cm)])
  }
  res <- lapply(1:10, densify_chr)
  G <- do.call(rbind, lapply(res, `[[`, "geno"))
  markers <- rbindlist(lapply(res, `[[`, "markers"))
  message(sprintf("  complete matrix: %d markers x %d lines; NA=%d", nrow(G), ncol(G), sum(is.na(G))))
  saveRDS(list(G = G, markers = markers), file.path(TDIR, "geno_gwas_interp.rds"))
  # recode 0/1/2 -> A/M/C (complete, no N)
  cn <- matrix("N", nrow(G), ncol(G))
  cn[G == 0] <- "A"
  cn[round(G) == 1] <- "M"
  cn[G == 2] <- "C"
  hd <- data.table(
    `rs#` = rownames(G), alleles = "A/C", chrom = as.integer(markers$chr),
    pos = as.integer(pos_by[rownames(G)]), strand = "+", `assembly#` = NA, center = NA,
    protLSID = NA, assayLSID = NA, panelLSID = NA, QCcode = NA
  )
  hmp <- cbind(hd, as.data.table(cn))
  setnames(hmp, c(names(hd), colnames(G)))
  setorder(hmp, chrom, pos)
  fwrite(hmp, HMP, sep = "\t", quote = FALSE, na = "N")
  message(sprintf("  wrote %s", basename(HMP)))

  # Q = 5 PCs (prcomp on complete G) + STAM phenotype (covariates)
  Xm <- t(G)
  Xm <- Xm[, (colMeans(Xm^2) - colMeans(Xm)^2) > 0, drop = FALSE]
  pc <- prcomp(Xm, center = TRUE, scale. = FALSE)$x[, 1:5]
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
  message(sprintf("  wrote %s (%d taxa)", basename(PHENO), nrow(pt)))
}

if (stage %in% c("kinship", "all")) {
  message("kinship (Centered_IBS on complete matrix) ...")
  r <- tass(sprintf('-importGuess "%s" -KinshipPlugin -method Centered_IBS -endPlugin -export "%s" -exportType SqrMatrix', HMP, KIN))
  message(sprintf("  kinship rc=%d, %.1f s", r$rc, r$sec))
  timings[["kinship"]] <- r$sec
}
if (stage %in% c("mlm", "all")) {
  message("MLM (P3D, Q+K) on complete matrix ...")
  r <- tass(sprintf('-fork1 -h "%s" -fork2 -r "%s" -fork3 -k "%s" -combine4 -input1 -input2 -intersect -combine5 -input4 -input3 -mlm -mlmVarCompEst P3D -mlmCompressionLevel None -export "%s" -runfork1 -runfork2 -runfork3', HMP, PHENO, KIN, MLMOUT))
  message(sprintf("  MLM rc=%d, %.1f s; files: %s", r$rc, r$sec, paste(list.files(TDIR, pattern = "mlm_interp"), collapse = ", ")))
  timings[["mlm"]] <- r$sec
}
if (length(timings)) {
  tt <- data.table(step = names(timings), sec = unlist(timings))
  fwrite(tt, file.path(OUT, "mlm_interp_timings.csv"))
  print(tt)
}
