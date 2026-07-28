#!/usr/bin/env Rscript
# Holland-style HMM grid search on the ZEAL population, through caller_grid().
# Grid = nir x germ x gert x p (Holland's emission levels) x rrate, where rrate is
# sampled ONE POINT PER ORDER OF MAGNITUDE from 1e-6 to 1e-2 = {1e-6,1e-5,1e-4,
# 1e-3,1e-2} (5 points), deliberately WIDER than Zhong/Holland's 3-point band
# (avg_r/2, avg_r, avg_r*2 ~ a 4x range around the map value; File_S16.py:168), so
# the full rrate response is mapped rather than just the neighbourhood of map r.
#
# Run on the ZEAL count-from-parents BC2S3 6-sib-pool sim (results/sim/zeal_pool),
# scored against the pooled-dose sim truth by the per-cell ANCESTRY-state mismatch
# (marker_dice, 1 - accuracy over all states). nnil is categorical: skim counts are
# hard-called once via the BC2S3 design-prior MAP (call_gt), then the whole grid is
# decoded on that fixed g. Chunked over the emission grid (score-and-discard, bounded
# memory), cores capped for <=16 GB, upfront ETA probe + per-batch ETA.
#
#   Rscript scripts/zeal_holland_grid.R
# Output: data/zeal/zeal_holland_grid.csv (per-config mismatch + fragment metrics)

suppressMessages({
  library(devtools)
  library(data.table)
  library(parallel)
})
devtools::load_all(path.expand("~/repos/nilhmm"), quiet = TRUE) # caller_grid, call_gt, breeding_prior
root <- here::here()
for (f in list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)) source(f) # marker_dice etc. (R/metrics.R)
source(file.path(root, "scripts/logging.R"))

MEM_CAP_GB <- 16L
N_CORES <- min(detectCores() - 2L, 8L)
N_CAL <- 120L # calibration line subset (full 49K grid; keep RAM bounded)
DESIGN <- "BC2S3"
SIM <- file.path(root, "results/sim/zeal_pool/zeal_pool_bc2s3_full.rds")

# emission grid (Holland levels) + wide rrate (5 points, one per order of magnitude)
EG <- as.data.table(expand.grid(
  nir = c(0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9),
  germ = c(1e-4, 1e-3, 1e-2), gert = c(1e-4, 1e-3, 1e-2),
  p = c(0.1, 0.25, 0.5, 0.75, 0.9)
))
RR <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
N_CFG <- nrow(EG) * length(RR)

sim <- readRDS(SIM)
design_prior <- breeding_prior(DESIGN)
ids <- head(sim$names, min(N_CAL, length(sim$names)))
idx <- match(ids, sim$names)
M <- nrow(sim$grid)
skim <- data.table(
  name = rep(ids, each = M), chr = rep(sim$grid$chr, length(idx)), pos = rep(sim$grid$pos, length(idx)),
  n_ref = as.integer(sim$n_ref[, idx]), n_alt = as.integer(sim$n_alt[, idx])
)
mr_hat <- skim[, mean(n_ref + n_alt == 0)] # realized skim missing rate -> fixed mr
EG[, mr := mr_hat]
skim[, g := {
  gg <- call_gt(n_ref, n_alt, prior = design_prior, error = 0.01)
  gg[is.na(gg)] <- 3L
  as.integer(gg)
}]
setorder(skim, name, chr, pos)
gdat <- skim[, .(name, chr, pos, g)]
truth <- sim$truth[name %in% ids]
grid_eval <- as.data.table(sim$grid)[, .(chr, pos)]
tr_sizes <- donor_block_sizes(truth)
log_info(
  "ZEAL Holland grid: %d emission x %d rrate = %d configs on %d cal lines x %d markers | mr=%.3f | cores=%d (<=%d GB)",
  nrow(EG), length(RR), N_CFG, length(ids), M, mr_hat, N_CORES, MEM_CAP_GB
)
log_info("rrate grid (1 point/order of magnitude): %s", paste(RR, collapse = ", "))

# score one batch of emission combos: caller_grid over the rrate grid, per-cell
# ancestry mismatch + fragment metrics vs the sim truth, then discard the segments.
score_batch <- function(eg_batch) {
  seg <- as.data.table(caller_grid(gdat,
    caller = "nnil", emission_grid = as.data.frame(eg_batch),
    rrate = RR, design = DESIGN, threads = N_CORES
  ))
  seg[, cfg := paste(nir, germ, gert, p, rrate, sep = "_")]
  cfgs <- unique(seg[, .(nir, germ, gert, p, r = rrate, cfg)])
  out <- rbindlist(lapply(seq_len(nrow(cfgs)), function(i) {
    called <- seg[cfg == cfgs$cfg[i]]
    mf <- marker_dice(called, truth, grid_eval)
    ff <- donor_fragment_dice(called, truth)
    data.table(
      nir = cfgs$nir[i], germ = cfgs$germ[i], gert = cfgs$gert[i], p = cfgs$p[i], r = cfgs$r[i],
      mismatch = 1 - mf$accuracy, donor_frag_dice = ff$dice,
      frag_ks = fragment_size_ks(donor_block_sizes(called), tr_sizes)
    )
  }))
  rm(seg)
  gc(FALSE)
  out
}

# upfront ETA probe: one emission combo x the rrate grid
tp <- Sys.time()
invisible(score_batch(EG[1]))
per_cfg <- as.numeric(difftime(Sys.time(), tp, units = "secs")) / length(RR)
log_info(
  "PROBE: %.2fs/config -> projected full grid ~%.1f min (%d configs, %d cores). Starting.",
  per_cfg, per_cfg * N_CFG / 60, N_CFG, N_CORES
)

# chunked run with per-batch ETA
BATCH <- 21L
batches <- split(seq_len(nrow(EG)), (seq_len(nrow(EG)) - 1L) %/% BATCH)
B <- length(batches)
t0 <- Sys.time()
res <- vector("list", B)
for (b in seq_len(B)) {
  res[[b]] <- score_batch(EG[batches[[b]]])
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  log_info(
    ">>> batch %d/%d done | elapsed %.1f min | avg %.2f min/batch | ETA ~%.1f min remaining",
    b, B, el, el / b, (el / b) * (B - b)
  )
}
score <- rbindlist(res)
fwrite(score, file.path(root, "data/zeal/zeal_holland_grid.csv"))
setorder(score, mismatch)
best <- score[1]
log_info("wrote data/zeal/zeal_holland_grid.csv (%d configs)", nrow(score))
log_info(
  "BEST mismatch=%.4f at nir=%g germ=%g gert=%g p=%g r=%.0e | fragDSC=%.3f KS=%.3f",
  best$mismatch, best$nir, best$germ, best$gert, best$p, best$r, best$donor_frag_dice, best$frag_ks
)
cat("\n=== mismatch vs rrate (marginal over emission), the wide-rrate curve ===\n")
print(score[, .(mismatch = round(mean(mismatch), 4)), by = r][order(r)])
