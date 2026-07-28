#!/usr/bin/env Rscript
# ZEAL BrB atlas NIR calibration at FIXED rigidity r=50 (the KS-optimal rigidity from
# fig_zeal_atlas_calibration.R). rrate/rigidity are held; the calibration knob is `nir`,
# the non-informative-marker rate of the gt emission (per the nnil clean plan: calibrate nir, not
# the duration). atlas = GOOGA-threshold categorical gt emission + rigidity duration, on the
# cassini pangene competitive counts. Criterion = introgression-size KS to the BC2S3 sim-latent
# fragment-size distribution rasterized onto the pangene grid (reuses the cached target).
#
# INCREMENTAL + RESUMABLE: each nir's segments cached, KS row appended immediately, re-run skips
# cached nir. NNIL_ATLASNIR_RECOMPUTE=1 clears caches.
#   Rscript scripts/fig_zeal_atlas_nir_calibration.R

suppressMessages({
  library(devtools)
  library(data.table)
  library(here)
  library(ggplot2)
  library(patchwork)
})
devtools::load_all(path.expand("~/repos/nilhmm"), quiet = TRUE)
source(here::here("R/metrics.R"))
source(here::here("scripts/logging.R"))

OUT <- here::here("results/sim/zeal_nil")
RIG_FIXED <- 50L
NIR_GRID <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9) # non-informative-marker rate (the knob)
GERM <- 0.05
GERT <- 0.10
P <- 0.5
MR <- 0.10 # other gt-emission params held at defaults
ATLAS_THRESH <- 0.95
ATLAS_HET <- 0.25
ATLAS_MIN_READS <- 5L
KEEP <- c("name", "chr", "start_bp", "end_bp", "state")
BASE <- 20
SEG_DIR <- file.path(OUT, "atlas_nir_seg_cache")
SWEEP_CSV <- file.path(OUT, "fig_zeal_atlas_nir_calibration_sweep.csv")
TARGET_RDS <- file.path(OUT, "atlas_target_sz.rds") # shared with the rigidity sweep
RECOMPUTE <- Sys.getenv("NNIL_ATLASNIR_RECOMPUTE") != ""
dir.create(SEG_DIR, showWarnings = FALSE, recursive = TRUE)
if (RECOMPUTE) {
  unlink(list.files(SEG_DIR, full.names = TRUE))
  unlink(SWEEP_CSV)
}
seg_cache <- function(v) file.path(SEG_DIR, sprintf("nir_%.3f.rds", v))

# ---- pangene grid + atlas input (same construction as the rigidity sweep) ----
gtp <- fread(here::here("data/ref/gene_to_pangene.tsv"))
gc <- fread(here::here("data/ref/b73_gene_coords.tsv"), header = FALSE, col.names = c("gene", "chr", "start", "end"))
gc[, chrn := as.integer(sub("chr", "", chr))]
pg <- merge(gtp[species == "B73", .(gene, pan_gene)], gc[!is.na(chrn), .(gene, chrn, start)], by = "gene")
pgchr <- pg[, .(nchr = uniqueN(chrn)), by = pan_gene]
pg <- pg[pan_gene %in% pgchr[nchr == 1L, pan_gene]]
pgpos <- pg[, .(chr = as.integer(chrn[1]), pos = as.integer(stats::median(start))), by = pan_gene]
pg_grid <- unique(pgpos[, .(chr, pos)])[order(chr, pos)]
pf <- list.files(here::here("data/skimsweep/brb/pangene"), pattern = "pangene_counts\\.tsv$", recursive = TRUE, full.names = TRUE)
brb <- rbindlist(lapply(pf, function(f) fread(f)[, .(name = sample, pan_gene, n_recur, n_donor)]))
brb <- merge(brb, pgpos, by = "pan_gene")
atlas_in <- brb[, .(name, chr, pos, n_ref = as.integer(n_recur), n_alt = as.integer(n_donor))]
setorder(atlas_in, name, chr, pos)
log_info(
  "[atlas-nir] %d BrB samples; fixed rigidity=%d; sweeping nir %s",
  uniqueN(atlas_in$name), RIG_FIXED, paste(NIR_GRID, collapse = ",")
)

# ---- sim-truth target (reuse the cached fragment sizes) ----------------------
if (file.exists(TARGET_RDS)) {
  target_sz <- readRDS(TARGET_RDS)
} else {
  mono <- file.path(OUT, "fig_zeal_atlas_calibration_cache.rds") # rigidity sweep stores target_sz here
  if (!file.exists(mono)) stop("[atlas-nir] target cache missing; run fig_zeal_atlas_calibration.R first")
  target_sz <- readRDS(mono)$target_sz
  saveRDS(target_sz, TARGET_RDS) # persist for reuse
}
log_info("[atlas-nir] target: %d BC2S3 sim donor fragments (cached)", length(target_sz))

# ---- incremental, resumable nir sweep at r=50 --------------------------------
threads <- min(parallel::detectCores() - 2L, 8L)
t0 <- Sys.time()
for (v in NIR_GRID) {
  f <- seg_cache(v)
  if (file.exists(f)) {
    log_info("[atlas-nir] nir=%.3f cached", v)
    s <- as.data.table(readRDS(f))
  } else {
    s <- as.data.table(call_ancestry(as.data.frame(atlas_in),
      caller = "atlas", rigidity = RIG_FIXED,
      nir = v, germ = GERM, gert = GERT, p = P, mr = MR,
      atlas_thresh = ATLAS_THRESH, atlas_het = ATLAS_HET, atlas_min_reads = ATLAS_MIN_READS,
      design = "BC2S3", threads = threads
    ))
    s <- s[, ..KEEP]
    s$nir <- v
    saveRDS(s, f)
    log_info("[atlas-nir] nir=%.3f decoded (%.0fs); cached", v, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }
  sz <- donor_block_sizes(s[, ..KEEP])
  row <- data.table(
    nir = v, ks = fragment_size_ks(sz, target_sz),
    n_donor_blocks = length(sz), median_block_mb = if (length(sz)) stats::median(sz) else NA_real_
  )
  fwrite(row, SWEEP_CSV, append = file.exists(SWEEP_CSV))
  log_info("[atlas-nir] nir=%.3f KS=%.3f appended", v, row$ks)
}

# ---- figure ------------------------------------------------------------------
sweep <- unique(fread(SWEEP_CSV), by = "nir", fromLast = TRUE)[order(nir)]
nir_star <- sweep$nir[which.min(sweep$ks)]
log_info("[atlas-nir] KS-optimal nir = %.3f (D=%.3f) at rigidity=%d", nir_star, min(sweep$ks, na.rm = TRUE), RIG_FIXED)
col_atlas <- "#009E73"
p_A <- ggplot(sweep, aes(nir, ks)) +
  geom_line(linewidth = 0.9, colour = col_atlas) +
  geom_point(size = 1.8, colour = col_atlas) +
  geom_point(data = sweep[which.min(ks)], size = 4, shape = 21, fill = "white", stroke = 1.2, colour = col_atlas) +
  expand_limits(y = 0) +
  labs(
    x = "atlas nir (non-informative-marker rate)", y = "introgression-size KS distance (D)",
    title = sprintf("BrB atlas nir calibration (r=%d fixed)\nKS-optimal nir=%.2f", RIG_FIXED, nir_star)
  ) +
  theme_bw(base_size = BASE) +
  theme(aspect.ratio = 1)
lab_t <- "BC2S3 sim latent (pangene grid)"
lab_star <- sprintf("atlas calls (nir=%.2f, KS-argmin)", nir_star)
lab_07 <- "atlas calls (nir=0.70, donor-calculated)"
atlas_star <- donor_block_sizes(as.data.table(readRDS(seg_cache(nir_star)))[, ..KEEP])
atlas_07 <- donor_block_sizes(as.data.table(readRDS(seg_cache(0.7)))[, ..KEEP])
ecdf_dt <- rbindlist(list(
  data.table(size_mb = target_sz, series = lab_t),
  data.table(size_mb = atlas_star, series = lab_star),
  data.table(size_mb = atlas_07, series = lab_07)
))
p_B <- ggplot(ecdf_dt, aes(size_mb, colour = series)) +
  stat_ecdf(linewidth = 1) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100)) +
  scale_colour_manual(values = setNames(c("black", col_atlas, "#D55E00"), c(lab_t, lab_star, lab_07)), name = NULL) +
  labs(x = "introgression size (Mb)", y = "ECDF", title = sprintf("atlas (real BrB, r=%d) vs sim target", RIG_FIXED)) +
  theme_bw(base_size = BASE) +
  theme(
    aspect.ratio = 1, legend.position = c(0.02, 0.99), legend.justification = c(0, 1),
    legend.background = element_rect(fill = "transparent", colour = NA)
  )
fig <- (p_A | p_B) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 22, face = "bold"))
ggsave(file.path(OUT, "fig_zeal_atlas_nir_calibration.png"), fig, width = 12, height = 6, dpi = 150)
log_info("[atlas-nir] wrote fig_zeal_atlas_nir_calibration.png")
