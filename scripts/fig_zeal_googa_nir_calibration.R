#!/usr/bin/env Rscript
# ZEAL BrB googa NIR calibration with rrate PINNED to the pangene-grid map-derived rate.
# googa = GOOGA-threshold categorical gt emission + GEOMETRIC duration (the faithful GOOGA/Veltsos
# caller; the geometric sibling of atlas). Per the established principle: rrate is physical, pinned
# to the map (interpolate the native TeoNAM v5 SNP50K cM onto the pangene positions, Holland
# 2L/(100M) form); the emission knob is `nir` (swept here to show it is a flat plateau, then pinned
# to the donor-calculated ~0.7, cf. the atlas nir figure). Criterion = introgression-size KS to the
# BC2S3 sim-latent fragment sizes on the pangene grid (reuses the cached target).
# INCREMENTAL + RESUMABLE. NNIL_GOOGANIR_RECOMPUTE=1 clears caches.
#   Rscript scripts/fig_zeal_googa_nir_calibration.R

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
NIR_GRID <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9)
GERM <- 0.05
GERT <- 0.10
P <- 0.5
MR <- 0.10
ATLAS_THRESH <- 0.95
ATLAS_HET <- 0.25
ATLAS_MIN_READS <- 5L
KEEP <- c("name", "chr", "start_bp", "end_bp", "state")
BASE <- 20
SEG_DIR <- file.path(OUT, "googa_nir_seg_cache")
SWEEP_CSV <- file.path(OUT, "fig_zeal_googa_nir_calibration_sweep.csv")
TARGET_RDS <- file.path(OUT, "atlas_target_sz.rds")
RECOMPUTE <- Sys.getenv("NNIL_GOOGANIR_RECOMPUTE") != ""
dir.create(SEG_DIR, showWarnings = FALSE, recursive = TRUE)
if (RECOMPUTE) {
  unlink(list.files(SEG_DIR, full.names = TRUE))
  unlink(SWEEP_CSV)
}
seg_cache <- function(v) file.path(SEG_DIR, sprintf("nir_%.3f.rds", v))

# ---- pangene grid + googa input (same construction as atlas) -----------------
gtp <- fread(here::here("data/ref/gene_to_pangene.tsv"))
gc <- fread(here::here("data/ref/b73_gene_coords.tsv"), header = FALSE, col.names = c("gene", "chr", "start", "end"))
gc[, chrn := as.integer(sub("chr", "", chr))]
pg <- merge(gtp[species == "B73", .(gene, pan_gene)], gc[!is.na(chrn), .(gene, chrn, start)], by = "gene")
pgchr <- pg[, .(nchr = uniqueN(chrn)), by = pan_gene]
pg <- pg[pan_gene %in% pgchr[nchr == 1L, pan_gene]]
pgpos <- pg[, .(chr = as.integer(chrn[1]), pos = as.integer(stats::median(start))), by = pan_gene]
pf <- list.files(here::here("data/skimsweep/brb/pangene"), pattern = "pangene_counts\\.tsv$", recursive = TRUE, full.names = TRUE)
brb <- rbindlist(lapply(pf, function(f) fread(f)[, .(name = sample, pan_gene, n_recur, n_donor)]))
brb <- merge(brb, pgpos, by = "pan_gene")
googa_in <- brb[, .(name, chr, pos, n_ref = as.integer(n_recur), n_alt = as.integer(n_donor))]
setorder(googa_in, name, chr, pos)

# ---- rrate PINNED to the pangene-grid map-derived rate -----------------------
# interpolate native TeoNAM v5 SNP50K cM onto the decoded pangene grid, Holland 2L/(100M).
cmmap <- fread(here::here("data/zeal/markers_snp50k_cm.tsv"))[, .(chr = as.integer(chr), pos = as.integer(pos), cm = as.numeric(cm))]
grid <- unique(googa_in[, .(chr, pos)])[order(chr, pos)]
grid[, cm := as.numeric(NA)]
for (cc in unique(grid$chr)) {
  m <- cmmap[chr == cc][order(pos)]
  i <- grid$chr == cc
  grid$cm[i] <- stats::approx(m$pos, m$cm, xout = grid$pos[i], rule = 2)$y
}
Lsum <- sum(vapply(split(grid$cm, grid$chr), function(x) max(x) - min(x), numeric(1)))
Mgrid <- nrow(grid)
MAP_RRATE <- 2 * Lsum / (100 * Mgrid)
log_info(
  "[googa-nir] pangene grid %d markers; total map %.0f cM; PINNED rrate=%.3e; sweeping nir %s",
  Mgrid, Lsum, MAP_RRATE, paste(NIR_GRID, collapse = ",")
)

# ---- sim-truth target (reuse cache; fall back to monolithic) -----------------
if (file.exists(TARGET_RDS)) {
  target_sz <- readRDS(TARGET_RDS)
} else {
  mono <- file.path(OUT, "fig_zeal_atlas_calibration_cache.rds")
  if (!file.exists(mono)) stop("[googa-nir] target cache missing; run fig_zeal_atlas_calibration.R first")
  target_sz <- readRDS(mono)$target_sz
  saveRDS(target_sz, TARGET_RDS)
}
log_info("[googa-nir] target: %d BC2S3 sim donor fragments (cached)", length(target_sz))

# ---- incremental, resumable nir sweep at pinned rrate ------------------------
threads <- min(parallel::detectCores() - 2L, 8L)
t0 <- Sys.time()
for (v in NIR_GRID) {
  f <- seg_cache(v)
  if (file.exists(f)) {
    log_info("[googa-nir] nir=%.3f cached", v)
    s <- as.data.table(readRDS(f))
  } else {
    s <- as.data.table(call_ancestry(as.data.frame(googa_in),
      caller = "googa", rrate = MAP_RRATE,
      nir = v, germ = GERM, gert = GERT, p = P, mr = MR,
      atlas_thresh = ATLAS_THRESH, atlas_het = ATLAS_HET, atlas_min_reads = ATLAS_MIN_READS,
      design = "BC2S3", threads = threads
    ))
    s <- s[, ..KEEP]
    s$nir <- v
    saveRDS(s, f)
    log_info("[googa-nir] nir=%.3f decoded (%.0fs); cached", v, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }
  sz <- donor_block_sizes(s[, ..KEEP])
  row <- data.table(
    nir = v, rrate = MAP_RRATE, ks = fragment_size_ks(sz, target_sz),
    n_donor_blocks = length(sz), median_block_mb = if (length(sz)) stats::median(sz) else NA_real_
  )
  fwrite(row, SWEEP_CSV, append = file.exists(SWEEP_CSV))
  log_info("[googa-nir] nir=%.3f KS=%.3f appended", v, row$ks)
}

# ---- figure ------------------------------------------------------------------
sweep <- unique(fread(SWEEP_CSV), by = "nir", fromLast = TRUE)[order(nir)]
nir_star <- sweep$nir[which.min(sweep$ks)]
log_info("[googa-nir] KS-optimal nir = %.3f (D=%.3f) at pinned rrate=%.3e", nir_star, min(sweep$ks, na.rm = TRUE), MAP_RRATE)
col_g <- "#0072B2"
p_A <- ggplot(sweep, aes(nir, ks)) +
  geom_line(linewidth = 0.9, colour = col_g) +
  geom_point(size = 1.8, colour = col_g) +
  geom_point(data = sweep[which.min(ks)], size = 4, shape = 21, fill = "white", stroke = 1.2, colour = col_g) +
  expand_limits(y = 0) +
  labs(
    x = "googa nir (non-informative-marker rate)", y = "introgression-size KS distance (D)",
    title = sprintf("BrB googa nir calibration\nrrate pinned to map (%.2e); KS-opt nir=%.2f", MAP_RRATE, nir_star)
  ) +
  theme_bw(base_size = BASE) +
  theme(aspect.ratio = 1)
lab_t <- "BC2S3 sim latent (pangene grid)"
lab_star <- sprintf("googa calls (nir=%.2f, KS-argmin)", nir_star)
lab_07 <- "googa calls (nir=0.70, donor-calculated)"
sz_star <- donor_block_sizes(as.data.table(readRDS(seg_cache(nir_star)))[, ..KEEP])
sz_07 <- donor_block_sizes(as.data.table(readRDS(seg_cache(0.7)))[, ..KEEP])
ecdf_dt <- rbindlist(list(
  data.table(size_mb = target_sz, series = lab_t),
  data.table(size_mb = sz_star, series = lab_star),
  data.table(size_mb = sz_07, series = lab_07)
))
p_B <- ggplot(ecdf_dt, aes(size_mb, colour = series)) +
  stat_ecdf(linewidth = 1) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100)) +
  scale_colour_manual(values = setNames(c("black", col_g, "#D55E00"), c(lab_t, lab_star, lab_07)), name = NULL) +
  labs(x = "introgression size (Mb)", y = "ECDF", title = "googa (real BrB, rrate=map) vs sim target") +
  theme_bw(base_size = BASE) +
  theme(
    aspect.ratio = 1, legend.position = c(0.02, 0.99), legend.justification = c(0, 1),
    legend.background = element_rect(fill = "transparent", colour = NA)
  )
fig <- (p_A | p_B) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 22, face = "bold"))
ggsave(file.path(OUT, "fig_zeal_googa_nir_calibration.png"), fig, width = 12, height = 6, dpi = 150)
log_info("[googa-nir] wrote fig_zeal_googa_nir_calibration.png")
