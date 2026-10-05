#!/usr/bin/env Rscript
# =============================================================================
# ZEAL PREDN_CLY — hyperspectral-predicted leaf nitrogen -> SpATS genotype BLUEs
# -----------------------------------------------------------------------------
# Predicted_N is NOT a measured trait: it is a per-scan prediction from
# hyperspectral reflectance, carrying its own model-space diagnostic
# (M_Distance = Mahalanobis distance of the spectrum from the calibration set).
# Scored on the CLY25 B5 block-2 grid only (921 of the field's 2,898 plots),
# 5 scans per plot, one plot per NIL.
#
# WHY THIS IS A SEPARATE SCRIPT (and not a trait added to zeal_spats_blues.R):
#   zeal_spats_blues.R reads trait columns straight out of the CLY25/CLY23
#   fieldbook workbooks for a hardcoded trait list, assuming one row per plot
#   with the trait ready to fit. Predicted_N arrives from a different source
#   (the BZeaPheno repo) at SCAN level and needs quality filtering and
#   plot-level aggregation before any spatial model can see it. That script is
#   also top-to-bottom executable, not a library, so sourcing it would
#   regenerate all 11 of its traits. The SpATS recipe below is deliberately the
#   same one (genotype FIXED -> BLUEs, PSANOVA(Row, Range) surface, row/col
#   factors random) so PREDN_CLY is comparable to its CLY siblings.
#
# PIPELINE
#   1. drop scans with M_Distance > M_MAX (default 3.5)
#   2. drop non-finite / non-positive predictions (a guard; none present as of
#      the 2025 data, where Predicted_N spans 2.136-4.931)
#   3. average surviving scans -> one value per plot
#   4. ONE SpATS fit over the block-2 grid (it is a single contiguous
#      rectangle, Row 18-45 x Range 2-35, 96.7% occupied; B73 appears in 125
#      plots spanning the full extent, which anchors the surface)
#   5. predict(which = "Genotype") -> BLUEs
#
# CHOOSING M_MAX: within-plot SD of Predicted_N rises with M_Distance but
# PLATEAUS above M = 3 (median within-plot SD 0.152 for M in 2-3, 0.151 for
# 3-5, 0.149 above 5). A cutoff below 3.5 therefore discards plots that are no
# less repeatable than ones it keeps, while costing whole genotypes (3.0 loses
# 15 genotypes, 2.0 loses 58). At 3.5: 141 scans (3.1%) and 5 genotypes drop.
# Note M_Distance is NOT independent of the trait (r = -0.24 with Predicted_N),
# so filtering trims both tails of the phenotype distribution — a reason to
# filter gently rather than strictly.
#
# In : data/zeal/B5_block2_predictedN.csv   (staged from the BZeaPheno repo,
#        data/B5_block2_predictedN.csv; cols Plot, Row, Range, Predicted_N,
#        M_Distance, Genotype, Accession, Taxa). Override with PREDN_RAW.
# Out: data/zeal/pheno_predn_cly_blue.csv   (Genotype, PREDN_cly25, PREDN_CLY_mean)
#      results/sim/zeal/zeal_predn_qc.csv   (per-plot scan counts + what was dropped)
# Env: M_MAX (3.5), PREDN_RAW
# Run: Rscript scripts/zeal_predn_pheno.R
# Next: TRAITS=PREDN_CLY NPERM=1000 NCORES=6 Rscript scripts/zeal_rqtl_scan.R
# =============================================================================
suppressMessages({
  library(here)
  library(data.table)
  library(SpATS)
})
source(here("scripts/logging.R"))

M_MAX <- as.numeric(Sys.getenv("M_MAX", "3.5")) # Mahalanobis model-space cutoff
RAW <- Sys.getenv("PREDN_RAW", here("data/zeal/B5_block2_predictedN.csv"))
ENV <- "CLY" # Clayton, the temperate macro-environment
FIELD <- "cly25" # the single field these scans come from
TAG <- paste0("PREDN_", ENV) # native-case trait token
OUT <- here(sprintf("data/zeal/pheno_%s_blue.csv", tolower(TAG)))
QC <- here("results/sim/zeal/zeal_predn_qc.csv")
dir.create(dirname(QC), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(RAW)) {
  stop(sprintf(
    "missing input: %s\nStage it from the BZeaPheno repo (data/B5_block2_predictedN.csv), or set PREDN_RAW.",
    RAW
  ))
}

# ---- 1-2. read + quality filters -------------------------------------------
d <- fread(RAW)
n0 <- nrow(d)
log_info("read %d scans | %d plots | %d genotypes", n0, uniqueN(d$Plot), uniqueN(d$Genotype))

bad <- d[!is.finite(Predicted_N) | Predicted_N <= 0, .N]
if (bad > 0) log_warn("dropping %d scans with non-finite / non-positive Predicted_N", bad)
d <- d[is.finite(Predicted_N) & Predicted_N > 0]

far <- d[M_Distance > M_MAX, .N]
dropped_geno <- setdiff(unique(d$Genotype), unique(d[M_Distance <= M_MAX]$Genotype))
d <- d[M_Distance <= M_MAX]
log_info(
  "M_Distance <= %.2f: dropped %d scans (%.2f%%) -> %d scans, %d plots, %d genotypes",
  M_MAX, far, 100 * far / n0, nrow(d), uniqueN(d$Plot), uniqueN(d$Genotype)
)
if (length(dropped_geno)) {
  log_warn("%d genotypes lost entirely (all scans above cutoff): %s",
    length(dropped_geno), paste(head(dropped_geno, 10), collapse = ", ")
  )
}

# ---- 3. scans -> plot means ------------------------------------------------
# One plot per NIL (B73 is the exception, 125 plots), so the plot mean IS the
# genotype's observation; the two-stage and one-stage averages agree to 4e-3.
plots <- d[, .(
  PREDN = mean(Predicted_N),
  n_scans = .N,
  sd_scans = ifelse(.N > 1, sd(Predicted_N), NA_real_),
  M_mean = mean(M_Distance)
), by = .(Plot, Genotype, Row, Range)]
thin <- plots[n_scans == 1, .N]
log_info("plot means: %d plots (%d rest on a single surviving scan)", nrow(plots), thin)

# ---- 4. SpATS over the single contiguous grid, genotype FIXED -> BLUEs ------
# Same model as zeal_spats_blues.R::fit_grid_trait(): the Rep term there is
# conditional on >1 level and this grid has one, so it is omitted here.
fit_dat <- data.frame(
  Genotype = factor(as.character(plots$Genotype)),
  Range = as.numeric(plots$Range),
  Row = as.numeric(plots$Row),
  y = as.numeric(plots$PREDN)
)
fit_dat$Rf <- factor(fit_dat$Range)
fit_dat$Cf <- factor(fit_dat$Row)
log_info(
  "SpATS: %d plots, %d genotypes, grid Row %d-%d x Range %d-%d",
  nrow(fit_dat), nlevels(fit_dat$Genotype),
  min(plots$Row), max(plots$Row), min(plots$Range), max(plots$Range)
)

m <- SpATS(
  response = "y", genotype = "Genotype", genotype.as.random = FALSE,
  spatial = ~ PSANOVA(Row, Range), random = ~ Rf + Cf,
  data = fit_dat, control = list(monitoring = 0, maxit = 50)
)

# ---- 5. genotype BLUEs -----------------------------------------------------
pr <- as.data.table(predict(m, which = "Genotype"))
vc <- grep("predicted", names(pr), ignore.case = TRUE, value = TRUE)[1]
bl <- pr[!is.na(get(vc)), .(Genotype = as.character(Genotype), blue = get(vc))]
setnames(bl, "blue", paste0("PREDN_", FIELD))
bl[, (paste0(TAG, "_mean")) := get(paste0("PREDN_", FIELD))]
setorder(bl, Genotype)
fwrite(bl, OUT)
log_info(
  "wrote %s | %d genotype BLUEs (mean %.3f, sd %.3f, range %.3f-%.3f)",
  basename(OUT), nrow(bl),
  mean(bl[[3]]), sd(bl[[3]]), min(bl[[3]]), max(bl[[3]])
)

# ---- QC companion ----------------------------------------------------------
qc <- plots[, .(Plot, Genotype, Row, Range, n_scans, sd_scans, M_mean, PREDN)]
setorder(qc, n_scans, -M_mean)
fwrite(qc, QC)
fwrite(
  data.table(Genotype = dropped_geno, reason = sprintf("all scans M_Distance > %.2f", M_MAX)),
  sub("[.]csv$", "_dropped_genotypes.csv", QC)
)
log_info("wrote QC: %s (+ _dropped_genotypes.csv)", basename(QC))
