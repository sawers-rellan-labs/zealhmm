#!/usr/bin/env Rscript
# ZEAL BrB atlas RIGIDITY calibration on the pangene grid -- the BrB-track analog of the Skim
# rtiger rigidity figure. atlas = GOOGA-threshold categorical gt emission + rigidity duration,
# run on the cassini pangene COMPETITIVE counts (n_recur = recurrent/B73, n_donor = donor).
# Dispatcher: RNA/ASE competitive counts -> categorical via GOOGA threshold -> ATLAS (rigid);
# never the BetaBinomial (red-x route). Rigidity is non-physical (min run in pangenes), so it is
# calibrated on the founder-sim fragment-length distribution:
#   target = BC2S3 sim latent ancestry (zeal_nil_bc2s3_full.rds $truth, native TeoNAM v5 map)
#            rasterized onto the pangene grid -> donor-fragment-size distribution.
#   sweep atlas rigidity on the REAL BrB pangene counts, minimize introgression-size KS to target.
#
# INCREMENTAL + RESUMABLE (atlas rigidity is O(T*(3r)^2) via state expansion, so high r is slow):
#   * each rigidity's segments are cached to atlas_seg_cache/rig_<r>.rds as it finishes;
#   * the KS row is APPENDED to the sweep CSV immediately (partial output always on disk);
#   * a re-run SKIPS any rigidity already cached. NNIL_ATLAS_RECOMPUTE=1 clears caches + CSV.
# Feasibility: atlas expands each macro-state into r sub-states, so a chromosome needs m >= 2r
#   markers; (sample,chr) chains failing 2r>m are dropped for that rigidity (logged).
#   Rscript scripts/fig_zeal_atlas_calibration.R

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
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
RIG_GRID <- as.integer(c(2, 3, 5, 7, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300, 500)) # atlas min-run (pangenes)
ATLAS_THRESH <- 0.95
ATLAS_HET <- 0.25
ATLAS_MIN_READS <- 5L
KEEP <- c("name", "chr", "start_bp", "end_bp", "state")
BASE <- 20
SEG_DIR <- file.path(OUT, "atlas_seg_cache")
SWEEP_CSV <- file.path(OUT, "fig_zeal_atlas_calibration_sweep.csv")
TARGET_RDS <- file.path(OUT, "atlas_target_sz.rds")
RECOMPUTE <- Sys.getenv("NNIL_ATLAS_RECOMPUTE") != ""
dir.create(SEG_DIR, showWarnings = FALSE)
if (RECOMPUTE) {
  unlink(list.files(SEG_DIR, full.names = TRUE))
  unlink(c(SWEEP_CSV, TARGET_RDS))
}
seg_cache <- function(r) file.path(SEG_DIR, sprintf("rig_%04d.rds", r))

# ---- 1. pangene grid: pan_gene -> B73 gene -> (chr, pos) ---------------------
gtp <- fread(here::here("data/ref/gene_to_pangene.tsv")) # gene, pan_gene, species
gc <- fread(here::here("data/ref/b73_gene_coords.tsv"),
  header = FALSE,
  col.names = c("gene", "chr", "start", "end")
)
gc[, chrn := as.integer(sub("chr", "", chr))]
pg <- merge(gtp[species == "B73", .(gene, pan_gene)], gc[!is.na(chrn), .(gene, chrn, start)], by = "gene")
pgchr <- pg[, .(nchr = uniqueN(chrn)), by = pan_gene] # drop paralog pan_genes spanning >1 chr
pg <- pg[pan_gene %in% pgchr[nchr == 1L, pan_gene]]
pgpos <- pg[, .(chr = as.integer(chrn[1]), pos = as.integer(stats::median(start))), by = pan_gene]
pg_grid <- unique(pgpos[, .(chr, pos)])[order(chr, pos)]
log_info("[atlas] pangene grid: %d anchored pan_genes on %d chr", nrow(pgpos), uniqueN(pgpos$chr))

# ---- 2. real BrB pangene competitive counts -> atlas input ------------------
pf <- list.files(here::here("data/skimsweep/brb/pangene"),
  pattern = "pangene_counts\\.tsv$", recursive = TRUE, full.names = TRUE
)
brb <- rbindlist(lapply(pf, function(f) fread(f)[, .(name = sample, pan_gene, n_recur, n_donor)]))
brb <- merge(brb, pgpos, by = "pan_gene")
atlas_in <- brb[, .(name, chr, pos, n_ref = as.integer(n_recur), n_alt = as.integer(n_donor))]
setorder(atlas_in, name, chr, pos)
nchain <- uniqueN(atlas_in[, .(name, chr)])
log_info(
  "[atlas] %d BrB samples x up to %d anchored pangene markers; %d (sample,chr) chains",
  uniqueN(atlas_in$name), uniqueN(atlas_in[, .(chr, pos)]), nchain
)

# ---- 3. sim-truth target on the pangene grid (cached) -----------------------
if (file.exists(TARGET_RDS)) {
  target_sz <- readRDS(TARGET_RDS)
  log_info("[atlas] target: reused cache (%d BC2S3 sim donor fragments)", length(target_sz))
} else {
  sim <- readRDS(file.path(OUT, "zeal_nil_bc2s3_full.rds"))
  N_SIM <- as.integer(Sys.getenv("NNIL_ATLAS_NSIM", "300"))
  set.seed(1)
  sub <- sim$names[sort(sample(seq_along(sim$names), min(N_SIM, length(sim$names))))]
  sim_truth <- as.data.table(sim$truth)[name %in% sub]
  tr <- rasterize_named(sim_truth, pg_grid)
  truth_long <- rbindlist(lapply(names(tr), function(nm) {
    data.table(name = nm, chr = pg_grid$chr, pos = pg_grid$pos, state = tr[[nm]])
  }))[!is.na(state)]
  truth_seg <- as.data.table(to_segments(as.data.frame(truth_long)))
  target_sz <- donor_block_sizes(truth_seg[, ..KEEP])
  saveRDS(target_sz, TARGET_RDS)
  log_info(
    "[atlas] target: %d BC2S3 sim donor fragments on the pangene grid (%d lines)",
    length(target_sz), length(sub)
  )
}

# ---- 4. incremental, resumable rigidity sweep -------------------------------
threads <- min(parallel::detectCores() - 2L, 8L)
t0 <- Sys.time()
for (r in RIG_GRID) {
  f <- seg_cache(r)
  if (file.exists(f)) {
    log_info("[atlas] rigidity=%d cached, skipping decode", r)
    s <- as.data.table(readRDS(f))
  } else {
    feas <- atlas_in[, if (.N >= 2L * r) .SD, by = .(name, chr)]
    ndrop <- nchain - uniqueN(feas[, .(name, chr)])
    if (!nrow(feas)) {
      log_info("[atlas] rigidity=%d: ALL %d chains infeasible (2r>m); skipped", r, nchain)
      next
    }
    s <- as.data.table(call_ancestry(as.data.frame(feas),
      caller = "atlas", rigidity = r,
      atlas_thresh = ATLAS_THRESH, atlas_het = ATLAS_HET, atlas_min_reads = ATLAS_MIN_READS,
      design = "BC2S3", threads = threads
    ))
    s <- s[, ..KEEP]
    s$rigidity <- r
    saveRDS(s, f)
    log_info(
      "[atlas] rigidity=%d decoded (%.0fs); %d/%d chains dropped (2r>m); cached",
      r, as.numeric(difftime(Sys.time(), t0, units = "secs")), ndrop, nchain
    )
  }
  sz <- donor_block_sizes(s[, ..KEEP])
  row <- data.table(
    rigidity = r, ks = fragment_size_ks(sz, target_sz),
    n_donor_blocks = length(sz), median_block_mb = if (length(sz)) stats::median(sz) else NA_real_
  )
  fwrite(row, SWEEP_CSV, append = file.exists(SWEEP_CSV)) # partial output persisted NOW
  log_info("[atlas] rigidity=%d KS=%.3f appended -> %s", r, row$ks, basename(SWEEP_CSV))
}

# ---- 5. figure from whatever is on disk (partial or full) -------------------
sweep <- unique(fread(SWEEP_CSV), by = "rigidity", fromLast = TRUE)[order(rigidity)]
rig_star <- sweep$rigidity[which.min(sweep$ks)]
log_info(
  "[atlas] KS-optimal rigidity so far = %d (D=%.3f) over %d rigidities on disk",
  rig_star, min(sweep$ks, na.rm = TRUE), nrow(sweep)
)

col_atlas <- "#009E73"
p_A <- ggplot(sweep, aes(rigidity, ks)) +
  geom_line(linewidth = 0.9, colour = col_atlas) +
  geom_point(size = 1.6, colour = col_atlas) +
  geom_point(data = sweep[which.min(ks)], size = 3.6, shape = 21, fill = "white", stroke = 1.1, colour = col_atlas) +
  scale_x_log10() +
  expand_limits(y = 0) +
  labs(
    x = "atlas rigidity (min run, pangenes)", y = "introgression-size KS distance (D)",
    title = sprintf("BrB atlas rigidity calibration\n(pangene grid; KS-opt so far r=%d)", rig_star)
  ) +
  theme_bw(base_size = BASE) +
  theme(aspect.ratio = 1)

atlas_star_sz <- donor_block_sizes(as.data.table(readRDS(seg_cache(rig_star)))[, ..KEEP])
ecdf_dt <- rbindlist(list(
  data.table(size_mb = target_sz, series = "BC2S3 sim latent\n(pangene grid)"),
  data.table(size_mb = atlas_star_sz, series = sprintf("atlas calls\n(r=%d)", rig_star))
))
p_B <- ggplot(ecdf_dt, aes(size_mb, colour = series)) +
  stat_ecdf(linewidth = 1) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100)) +
  scale_colour_manual(values = setNames(
    c("black", col_atlas),
    c("BC2S3 sim latent\n(pangene grid)", sprintf("atlas calls\n(r=%d)", rig_star))
  ), name = NULL) +
  labs(x = "introgression size (Mb)", y = "ECDF", title = "atlas (real BrB) vs BC2S3 sim target") +
  theme_bw(base_size = BASE) +
  theme(
    aspect.ratio = 1, legend.position = c(0.02, 0.99), legend.justification = c(0, 1),
    legend.background = element_rect(fill = "transparent", colour = NA)
  )

fig <- (p_A | p_B) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 22, face = "bold"))
ggsave(file.path(OUT, "fig_zeal_atlas_calibration.png"), fig, width = 12, height = 6, dpi = 150)
log_info("[atlas] wrote fig_zeal_atlas_calibration.png (%d rigidities)", nrow(sweep))
