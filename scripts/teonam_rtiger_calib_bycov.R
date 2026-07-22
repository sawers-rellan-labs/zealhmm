#!/usr/bin/env Rscript
# =============================================================================
# Per-coverage rtiger rigidity calibration on SIMULATED ground truth (118K grid).
#
# Why per-coverage: a fair coverage-degradation sweep needs each caller at its
# per-coverage optimum, not one rigidity tuned at a single depth. Ground truth is a
# SIMULATION (simulate_nil, known crossovers) -- NOT a real family's FSFHap calls
# (that would score a caller against another caller). The sim uses the SAME 118K v5
# marker grid as the experiment (rigidity is a marker count, so density must match)
# and the SAME read model as the sweep (.draw_counts, error=0.01).
#
# Method (per coverage lambda, reusing R/calibrate.R):
#   coarse sweep on a mean-density-centered powers-of-2 grid {8..512}  (markers/Mb
#   mean ~57 -> center 2^6=64), feasibility-filtered -> bracket_from_sweep ->
#   golden_refine, objective = donor-fragment Dice. Ties break toward the LARGER
#   rigidity ("err long": with linkage, over-merging costs less than over-fragmenting).
#   min_reads = 1 (decode only read-covered markers; ~2.6x faster at low coverage).
#
# Output (ephemeral sim, only small CSVs written):
#   results/sim/teonam/rtiger_calib_bycov.csv        (rigidity*(lambda) + Dice)
#   results/sim/teonam/rtiger_calib_bycov_sweep.csv  (full coarse Dice curves)
# Run: Rscript scripts/teonam_rtiger_calib_bycov.R [--smoke]   (smoke = lambda=1 only)
# =============================================================================
suppressMessages({
  library(data.table)
  devtools::load_all("/Users/fvrodriguez/repos/nilhmm", quiet = TRUE)
})
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)
for (f in list.files(file.path(ROOT, "R"), "\\.R$", full.names = TRUE)) source(f) # calibrate.R + simulate.R (.draw_counts)
OUT <- file.path(ROOT, "results/sim/teonam")
t0 <- Sys.time()
source(file.path(ROOT, "scripts/logging.R"))

SMOKE <- "--smoke" %in% commandArgs(TRUE)
SEED <- 12345L
N_RIL <- 200L
DESIGN <- "BC1S4"
READ_PARS <- list(pi_floor = 0, k_decay = 1, error = 0.01)
LAMBDAS <- if (SMOKE) c(1) else c(0.1, 0.2, 0.5, 1, 5, 10, 20, Inf)
THREADS <- max(1L, parallel::detectCores() - 2L)

# --- simulate BC1S4 ground truth on the experimental 118K v5 cM grid ----------
mc <- fread("data/teonam/markers_v5_gwas118k_cm.tsv")
setnames(mc, "pos_v5", "pos")
map <- data.frame(chr = as.integer(mc$chr), cm = as.numeric(mc$cm), bp = as.integer(mc$pos))
log_info("simulating %s n=%d on the %d-marker 118K grid (seed=%d) ...", DESIGN, N_RIL, nrow(map), SEED)
ts <- Sys.time()
truth_m <- as.data.table(simulate_nil(DESIGN, n = N_RIL, map = map, n_markers = nrow(map), seed = SEED))
truth <- as.data.table(to_segments(truth_m)) # common-schema segments for the metrics
grid <- unique(truth_m[, .(chr = as.integer(chr), pos = as.integer(pos))])
setorder(grid, chr, pos)
setkey(truth_m, name, chr, pos) # per-marker state lookup for the read draws
log_info(
  "  truth: %d segments, %d markers, %d lines; %.1f breakpoints/line (%.1f s)",
  nrow(truth), nrow(grid), uniqueN(truth_m$name), breakpoint_count(truth) / uniqueN(truth_m$name),
  as.numeric(difftime(Sys.time(), ts, units = "secs"))
)

# --- mean-density-centered powers-of-2 rigidity grid --------------------------
mpm <- mc[, .N, by = .(chr, mb = floor(pos / 1e6))]
center <- 2^round(log2(mean(mpm$N))) # mean markers/Mb ~57 -> 2^6 = 64
RIG_VALUES <- as.integer(center * 2^(-4:3)) # {4,8,16,32,64,128,256,512} (extra low octave)
log_info(
  "markers/Mb mean=%.0f -> center=%d; rigidity grid = %s",
  mean(mpm$N), center, paste(RIG_VALUES, collapse = ", ")
)
log_info(
  "=== rtiger per-coverage rigidity calibration: %d coverages x %d rigidities, %d threads ===",
  length(LAMBDAS), length(RIG_VALUES), THREADS
)

# --- read draws for one coverage (matches the sweep's .draw_counts model) -----
make_counts <- function(lambda) {
  set.seed(SEED + 1000L * which(LAMBDAS == lambda)) # per-coverage reproducible draw
  st <- truth_m$state
  if (is.infinite(lambda)) { # perfect-coverage ceiling: deterministic decisive counts
    p_alt <- c(0, 0.5, 1)[st + 1L]
    p_eff <- p_alt * (1 - READ_PARS$error) + (1 - p_alt) * READ_PARS$error
    n_alt <- as.integer(round(100L * p_eff))
    n_ref <- 100L - n_alt
  } else {
    ac <- .draw_counts(st, lambda, READ_PARS$pi_floor, READ_PARS$k_decay, READ_PARS$error)
    n_ref <- ac$ref
    n_alt <- ac$alt
  }
  data.table(
    name = truth_m$name, chr = as.integer(truth_m$chr), pos = as.integer(truth_m$pos),
    n_ref = n_ref, n_alt = n_alt
  )
}

# --- per-coverage calibration -------------------------------------------------
sweeps <- list()
best <- list()
t_loop <- Sys.time()
for (i in seq_along(LAMBDAS)) {
  lam <- LAMBDAS[i]
  covlab <- if (is.infinite(lam)) "Inf" else as.character(lam)
  d <- make_counts(lam)
  cov_frac <- mean(d$n_ref + d$n_alt > 0)
  vals <- feasible_rigidity(d, RIG_VALUES) # 2*r < min covered markers/chr
  t0 <- Sys.time()
  # ONE shared rtiger EM fit + a decode per rigidity (EM is rigidity-independent),
  # so no golden-refine (which would re-fit per probe). The integer grid + flat-top
  # Dice make sub-octave refinement pointless.
  sw <- sweep_calibrate(d, truth, grid,
    caller = "rtiger", values = vals,
    threads = THREADS, min_reads = 1L
  )
  el <- as.numeric(Sys.time() - t0, units = "mins")
  sw[, coverage := covlab]
  sweeps[[covlab]] <- sw
  # objective = donor-fragment FDR (MINIMIZE): Dice is monotone in rigidity here
  # (dominated by saturated recall) so it never brackets an optimum; FDR has a real
  # minimum (few spurious short fragments at low r, merged over-extended fragments at
  # high r). err-long: largest rigidity whose FDR is within TOL of the min (linkage ->
  # over-merging costs less than over-fragmenting; break near-ties toward longer).
  TOL <- 0.003
  fmin <- min(sw$donor_frag_FDR)
  rig_star <- sw[donor_frag_FDR <= fmin + TOL, max(value)]
  rig_argmin <- sw[donor_frag_FDR == fmin, value[1]]
  rc <- sw[value == rig_star, donor_marker_recall]
  best[[covlab]] <- data.table(
    coverage = covlab, cov_frac = round(cov_frac, 3),
    rigidity = as.integer(rig_star), rigidity_argmin = as.integer(rig_argmin),
    fdr = round(fmin, 4), recall = round(rc, 4),
    dice = round(sw[value == rig_star, donor_frag_dice], 4), mins = round(el, 1)
  )
  log_info(
    "  [%d/%d] lambda=%-4s (%.0f%% cov): rigidity* = %d (argmin-FDR %d), FDR=%.3f recall=%.3f (%.1f min)",
    i, length(LAMBDAS), covlab, 100 * cov_frac, as.integer(rig_star), as.integer(rig_argmin), fmin, rc, el
  )
  eltot <- as.numeric(difftime(Sys.time(), t_loop, units = "mins"))
  avg <- eltot / i
  log_info(
    ">>> %d/%d coverages done | elapsed %.1f min | avg %.1f min/cov | ETA ~%.1f min remaining",
    i, length(LAMBDAS), eltot, avg, avg * (length(LAMBDAS) - i)
  )
}

B <- rbindlist(best)
S <- rbindlist(sweeps)
if (!SMOKE) {
  fwrite(B, file.path(OUT, "rtiger_calib_bycov.csv"))
  fwrite(S, file.path(OUT, "rtiger_calib_bycov_sweep.csv"))
  log_info("wrote rtiger_calib_bycov.csv + _sweep.csv")
}
log_info(
  "=== DONE: %d coverages in %.1f min total ===", length(LAMBDAS),
  as.numeric(difftime(Sys.time(), t0, units = "mins"))
)
print(B)
