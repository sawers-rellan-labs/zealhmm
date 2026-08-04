#!/usr/bin/env Rscript
# Fragment-size cM-space validation of the five paint callers (linear-scale ECDF), the current-
# data analog of the legacy zealtiger fragment_size_cm_validation Part 3. Each caller runs at its
# recorded operating point (reference/zeal_skim_calibration_operating_points.csv) on the vary_skim
# coverage-sweep NILs (skim) / cassini pangene (BrB); its called donor blocks are converted Mb->cM
# via nilHMM::bp_to_cm (Hyman monotone Marey spline) fed the NATIVE TeoNAM v5 map
# (markers_snp50k_cm.tsv, 1559 cM -- NOT the bundled Ed Coe map). Reference = the BC2S3 sim latent
# donor-tract cM law (+ gamma fit). Linear cM x so the ~exponential BC2S3 law reads straight.
# Skim caller segments cached in paint_seg_cache/ (resumable); NNIL_CMVAL_RECOMPUTE=1 clears them.
#   Rscript scripts/fig_zeal_fragment_size_cm_validation.R

suppressMessages({
  library(devtools)
  library(data.table)
  library(here)
  library(ggplot2)
})
devtools::load_all(path.expand("~/repos/nilhmm"), quiet = TRUE)
source(here::here("R/metrics.R"))
source(here::here("scripts/logging.R"))
OUT <- here::here("results/sim/zeal_nil")
BASE <- 20
KEEP <- c("name", "chr", "start_bp", "end_bp", "state")
SEG_DIR <- file.path(OUT, "paint_seg_cache")
dir.create(SEG_DIR, showWarnings = FALSE)
# only clears binhmm.rds (the one still fit here); nnil/bbnil/rtiger are owned by build_paint_seg_cache.R
if (Sys.getenv("NNIL_CMVAL_RECOMPUTE") != "") unlink(file.path(SEG_DIR, "binhmm.rds"))
threads <- min(parallel::detectCores() - 2L, 8L)

# ---- native TeoNAM v5 Marey spline (bp -> cM), reused from nilHMM -------------
nmap <- fread(here::here("data/zeal/markers_snp50k_cm.tsv"))[, .(chr = as.integer(chr), bp = as.integer(pos), cm = as.numeric(cm))]
to_cm <- bp_to_cm(as.data.frame(nmap)) # Hyman monotone spline per chr
db_cm <- function(seg) { # donor-block cM lengths (>0)
  db <- .donor_blocks(as.data.table(seg)[, ..KEEP])
  if (!nrow(db)) {
    return(numeric(0))
  }
  cm <- to_cm(db$chr, db$end_bp) - to_cm(db$chr, db$start_bp)
  cm[is.finite(cm) & cm > 0]
}

# ---- skim binhmm on 1 Mb bins (the only caller still fit here; nnil/bbnil/rtiger are pooled in
#      scripts/build_paint_seg_cache.R). load_skim was retired with the fit-on-11 caller production.
load_bins <- function() {
  fs <- list.files(here::here("data/skimsweep/skim/bins"), pattern = "\\.tsv$", full.names = TRUE)
  rbindlist(lapply(fs, function(f) {
    b <- fread(f)
    data.table(
      name = b$SAMPLE, chr = as.integer(sub("chr", "", b$CONTIG)), pos = as.integer(b$BIN_POS),
      alt_freq = as.numeric(b$ALT_FREQ), start_bp = as.integer(b$BIN_START),
      end_bp = as.integer(b$BIN_END), ninf = as.integer(b$INFORMATIVE_VARIANT_COUNT)
    )
  }))
}
run_cached <- function(tag, fn) {
  f <- file.path(SEG_DIR, paste0(tag, ".rds"))
  if (file.exists(f)) {
    log_info("[cmval] %s cached", tag)
    return(as.data.table(readRDS(f)))
  }
  t0 <- Sys.time()
  s <- as.data.table(fn())[, ..KEEP]
  saveRDS(s, f)
  log_info("[cmval] %s decoded (%.0fs)", tag, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  s
}
# nnil/bbnil/rtiger segments are produced by scripts/build_paint_seg_cache.R as a SINGLE-POOL fit on
# the 330 dnarna cohort at each operating point (NOT fit here on the 11 coverage-sweep samples, which
# was an unrepresentative basis). This script now CONSUMES that cache; run build_paint_seg_cache.R
# first. binhmm stays fit here (1 Mb bins, a separate modality/input). See [[50k-set-terminology]].
require_cached <- function(tag) {
  f <- file.path(SEG_DIR, paste0(tag, ".rds"))
  if (!file.exists(f)) {
    stop(sprintf("[cmval] %s.rds missing -- run scripts/build_paint_seg_cache.R first (pool-330 fit)", tag))
  }
  as.data.table(readRDS(f))
}
nnil_seg <- require_cached("nnil")
bbnil_seg <- require_cached("bbnil")
rtiger_seg <- require_cached("rtiger")
binhmm_seg <- run_cached("binhmm", function() call_ancestry(as.data.frame(load_bins()), caller = "binhmm", design = "BC2S3"))
# ---- BrB callers from the nir-calibration caches (operating points) ----------
googa_seg <- as.data.table(readRDS(file.path(OUT, "googa_nir_seg_cache", "nir_0.200.rds")))
atlas_seg <- as.data.table(readRDS(file.path(OUT, "atlas_nir_seg_cache", "nir_0.100.rds")))

# ---- reference: BC2S3 sim latent donor-tract cM law + gamma fit --------------
sim_truth <- as.data.table(readRDS(file.path(OUT, "zeal_nil_bc2s3_full.rds"))$truth)
ref_cm <- db_cm(sim_truth)
gfit <- MASS::fitdistr(ref_cm / 100, "gamma")
k <- gfit$estimate["shape"]
lam <- gfit$estimate["rate"]
log_info(
  "[cmval] reference BC2S3 cM: n=%d mean=%.1f cM | gamma k=%.2f lambda=%.2f/M (mean %.1f cM)",
  length(ref_cm), mean(ref_cm), k, lam, 100 * k / lam
)

# ---- assemble + linear-cM ECDF -----------------------------------------------
L_ref <- "BC2S3 sim latent (truth)"
L <- c(
  nnil = "Skim-nnil (hard-call BC2S3, rrate=map, nir=0.70)",
  bbnil = "Skim-bbnil (rrate=map, fit_means)", rtiger = "Skim-rtiger (r=5)",
  binhmm = "Skim-binhmm (stay=0.995)", googa = "BrB-googa (rrate=map, nir=0.20)",
  atlas = "BrB-atlas (r=50, nir=0.10)"
)
dt <- rbindlist(list(
  data.table(cm = ref_cm, series = L_ref),
  data.table(cm = db_cm(nnil_seg), series = L[["nnil"]]),
  data.table(cm = db_cm(bbnil_seg), series = L[["bbnil"]]),
  data.table(cm = db_cm(rtiger_seg), series = L[["rtiger"]]),
  data.table(cm = db_cm(binhmm_seg), series = L[["binhmm"]]),
  data.table(cm = db_cm(googa_seg), series = L[["googa"]]),
  data.table(cm = db_cm(atlas_seg), series = L[["atlas"]])
))
lvl <- c(L_ref, L[["nnil"]], L[["bbnil"]], L[["rtiger"]], L[["binhmm"]], L[["googa"]], L[["atlas"]])
dt[, series := factor(series, levels = lvl)]
pal <- setNames(c("black", "#009E73", "#0072B2", "#56B4E9", "#999999", "#D55E00", "#E69F00"), lvl)
lty <- setNames(c("dotted", rep("solid", 6)), lvl)
ks <- dt[series != L_ref, .(D = fragment_size_ks(cm, ref_cm), n = .N, med = median(cm)), by = series]
log_info("[cmval] per-caller cM KS to BC2S3 truth + medians:")
print(ks, digits = 3)
xmax <- as.numeric(quantile(dt$cm, 0.99))

gcurve <- data.table(cm = seq(0, xmax, length.out = 300))
gcurve[, ecdf := pgamma(cm / 100, shape = k, rate = lam)]
p <- ggplot(dt, aes(cm, colour = series, linetype = series)) +
  stat_ecdf(linewidth = 1) +
  geom_line(data = gcurve, aes(cm, ecdf), inherit.aes = FALSE, colour = "grey55", linetype = "22", linewidth = 0.7) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_linetype_manual(values = lty, name = NULL) +
  coord_cartesian(xlim = c(0, xmax)) +
  labs(
    x = "donor introgression size (cM, native TeoNAM v5)", y = "ECDF",
    title = "Fragment-size cM validation of the six paint callers",
    subtitle = sprintf("linear cM  |  dotted black = BC2S3 sim latent truth\ngrey dashed = BC2S3 gamma fit (k=%.2f, mean %.1f cM)", k, 100 * k / lam)
  ) +
  theme_bw(base_size = BASE) +
  theme(
    legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
    legend.text = element_text(size = 11), legend.background = element_rect(fill = "transparent", colour = NA)
  )
ggsave(file.path(OUT, "fig_zeal_fragment_size_cm_validation.png"), p, width = 10, height = 8, dpi = 150)
fwrite(ks, file.path(OUT, "fig_zeal_fragment_size_cm_validation_ks.csv"))
log_info("[cmval] wrote fig_zeal_fragment_size_cm_validation.png")
