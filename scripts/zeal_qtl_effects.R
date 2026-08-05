#!/usr/bin/env Rscript
# scripts/zeal_qtl_effects.R
#
# DRAFT. Fills the two data gaps needed by the ZEAL 4-panel QTL summary figure:
#   1. per-QTL percent variance explained + additive effect (teosinte-B vs B73-A)
#   2. per-trait total QTL-model variance explained (full-model R2)
#
# For each trait: load the R/qtl bcsft cross and its confidence-interval peak
# table (zeal_<trait>_peaks_ci_taxon.csv, one peak per chromosome from the
# taxon-covariate joint scan), build a multi-QTL model at the peak markers
# (makeqtl on the peak cM positions), then fitqtl with a formula summing the
# QTL. From fitqtl we read the drop-one percent variance per QTL, the additive
# effect per QTL, and the full-model R2.
#
# Outputs (tidy CSVs, one dir results/sim/zeal/rqtl):
#   zeal_<trait>_qtl_effects.csv    trait,chr,pos_mb,marker,lod,pct_var,add_effect
#   zeal_qtl_variance_partition.csv trait,n_qtl,pct_var_qtl_model   (appended per trait)
#
# Robust/resumable: a trait whose fitqtl fails is logged and skipped; the run
# continues. Genotype codes in the cross are 1/2/3 = A/H/B where A = B73 recurrent,
# H = het, B = homozygous teosinte donor.
#
# NOTE: no SNP-heritability / taxon variance component is computed here; that is
# not readily available from fitqtl alone and is deliberately left out (not
# fabricated). Panel D therefore shows QTL-model R2 only.

suppressMessages(library(qtl))

# ---- config -----------------------------------------------------------------
TRAITS <- c("dta", "dts", "ph", "prolif") # 4 TeoNAM-overlap traits

# map lowercase trait key -> phenotype column name in the cross object
PHENO_COL <- c(dta = "DTA", dts = "DTS", ph = "PH", prolif = "Prolif")

RQTL_DIR <- "results/sim/zeal/rqtl"
PARTITION_CSV <- file.path(RQTL_DIR, "zeal_qtl_variance_partition.csv")

# start a fresh partition table for this run
if (file.exists(PARTITION_CSV)) invisible(file.remove(PARTITION_CSV))

ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
logmsg <- function(...) cat(sprintf("[%s] ", ts()), ..., "\n", sep = "")

# peak bp (v5) is parseable from the marker name "S<chr>_<bp>"
marker_bp <- function(m) as.numeric(sub("^S[0-9]+_", "", m))

# ---- per-trait effect estimation --------------------------------------------
fit_one_trait <- function(trait) {
  pcol <- PHENO_COL[[trait]]
  cross_path <- file.path(RQTL_DIR, sprintf("zeal_%s_cross.rds", trait))
  peaks_path <- file.path(RQTL_DIR, sprintf("zeal_%s_peaks_ci_taxon.csv", trait))

  logmsg("trait ", trait, " (pheno ", pcol, "): loading cross + peaks")
  cross <- readRDS(cross_path)
  peaks <- read.csv(peaks_path, stringsAsFactors = FALSE)
  if (nrow(peaks) == 0) stop("no peaks for ", trait)

  # genoprob for Haley-Knott fit (error.prob small, dense step not needed since
  # we evaluate at real peak markers via makeqtl)
  cross <- calc.genoprob(cross, step = 0, error.prob = 1e-4)

  # multi-QTL model at the ci peak cM positions (one per chromosome)
  qtl <- makeqtl(cross, chr = peaks$chr, pos = peaks$pos, what = "prob")
  form <- as.formula(paste("y ~", paste0("Q", seq_len(nrow(peaks)), collapse = " + ")))

  fq <- fitqtl(cross,
    pheno.col = pcol, qtl = qtl, formula = form,
    method = "hk", get.ests = TRUE, dropone = TRUE
  )
  sfq <- summary(fq)

  # full-model R2 (percent variance explained by the joint QTL model)
  full_pctvar <- sfq$result.full["Model", "%var"]

  # per-QTL drop-one percent variance. With >1 QTL fitqtl returns a drop-one
  # table with one row per QTL; with a single QTL there is no drop-one table so
  # we fall back to the full-model %var.
  if (nrow(peaks) > 1 && !is.null(sfq$result.drop)) {
    drop_tab <- sfq$result.drop
    pct_var <- drop_tab[, "%var"]
    # align drop-one rows to peaks by order (fitqtl preserves QTL order)
    pct_var <- pct_var[seq_len(nrow(peaks))]
  } else {
    pct_var <- rep(full_pctvar, nrow(peaks))
  }

  # additive effect per QTL: get.ests names the additive coefficient
  # "<qtl name>a" (and dominance "<qtl name>d"), where <qtl name> is the QTL
  # label from the peaks table (e.g. "9@6.5a"). a > 0 means the teosinte-B
  # homozygote raises the trait relative to the B73-A homozygote.
  ests <- fq$ests$ests
  add_effect <- vapply(seq_len(nrow(peaks)), function(i) {
    nm <- paste0(peaks$name[i], "a")
    if (nm %in% names(ests)) unname(ests[nm]) else NA_real_
  }, numeric(1))

  effects <- data.frame(
    trait = trait,
    chr = peaks$chr,
    pos_mb = marker_bp(peaks$marker) / 1e6,
    marker = peaks$marker,
    lod = peaks$lod,
    pct_var = round(as.numeric(pct_var), 4),
    add_effect = round(add_effect, 5),
    stringsAsFactors = FALSE
  )

  partition <- data.frame(
    trait = trait,
    n_qtl = nrow(peaks),
    pct_var_qtl_model = round(as.numeric(full_pctvar), 4),
    stringsAsFactors = FALSE
  )

  list(effects = effects, partition = partition)
}

# ---- run --------------------------------------------------------------------
logmsg("ZEAL QTL effects: traits = ", paste(TRAITS, collapse = ", "))

ok <- character(0)
failed <- character(0)

for (trait in TRAITS) {
  res <- tryCatch(fit_one_trait(trait), error = function(e) {
    logmsg("FAILED trait ", trait, ": ", conditionMessage(e))
    NULL
  })
  if (is.null(res)) {
    failed <- c(failed, trait)
    next
  }

  eff_path <- file.path(RQTL_DIR, sprintf("zeal_%s_qtl_effects.csv", trait))
  write.csv(res$effects, eff_path, row.names = FALSE)
  logmsg("wrote ", eff_path, " (", nrow(res$effects), " QTL)")

  # append the partition row (build-up across traits)
  write.table(res$partition, PARTITION_CSV,
    sep = ",", row.names = FALSE,
    col.names = !file.exists(PARTITION_CSV), append = file.exists(PARTITION_CSV)
  )
  ok <- c(ok, trait)
}

logmsg("done. succeeded: ", paste(ok, collapse = ", "))
if (length(failed)) logmsg("failed: ", paste(failed, collapse = ", "))
