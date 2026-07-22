#!/usr/bin/env Rscript
# =============================================================================
# Per-coverage nnil rrate calibration on SIMULATED ground truth (118K grid).
#
# The nnil analog of scripts/teonam_rtiger_calib_bycov.R. The tuned knob is `rrate`
# (per-marker recombination/transition rate of the geometric-duration HMM): low rrate
# -> long segments, high rrate -> fragmented. Unlike rtiger's integer rigidity, rrate
# is CONTINUOUS and nnil has no EM (emission fixed from `err` + the BC1S4 priors
# f_1/f_2), so we keep the golden-ratio refine (cheap here -- each probe is a decode,
# not a refit) and center the coarse grid on the density default
# rrate0 = 2*total_cM/(100*n_markers).
#
# Ground truth is a SIMULATION (simulate_nil, known crossovers) on the SAME 118K grid
# and read model (.draw_counts, error=0.01) as the sweep; objective = MIN donor-fragment
# FDR; min_reads = 1.
#
# Output (ephemeral sim, only small CSVs written):
#   results/sim/teonam/nnil_calib_bycov.csv        (rrate*(lambda) + FDR/recall)
#   results/sim/teonam/nnil_calib_bycov_sweep.csv  (full coarse FDR curves)
# Run: Rscript scripts/teonam_nnil_calib_bycov.R [--smoke]   (smoke = lambda=1 only)
# =============================================================================
suppressMessages({
  library(data.table)
  devtools::load_all("/Users/fvrodriguez/repos/nilhmm", quiet = TRUE)
})
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)
for (f in list.files(file.path(ROOT, "R"), "\\.R$", full.names = TRUE)) source(f) # calibrate.R (+ .draw_counts, log_grid)
OUT <- file.path(ROOT, "results/sim/teonam")
source(file.path(ROOT, "scripts/logging.R"))

SMOKE <- "--smoke" %in% commandArgs(TRUE)
SEED <- 12345L
N_RIL <- 200L
DESIGN <- "BC1S4"
READ_PARS <- list(pi_floor = 0, k_decay = 1, error = 0.01)
LAMBDAS <- if (SMOKE) c(1) else c(0.1, 0.2, 0.5, 1, 5, 10, 20, Inf)
# caller_sweep now segments in-worker (nilhmm), so the master no longer rbinds the
# full marker matrix -- decode workers hold one sample transiently and the master
# aggregates only compact segments. Safe to use the rtiger thread convention.
THREADS <- as.integer(Sys.getenv("NNIL_THREADS", as.character(max(1L, parallel::detectCores() - 2L))))

# nnil start priors for BC1S4 from breeding_prior() -- same as 02_calibrate.R
EXP <- breeding_prior("BC1S4") # (n_bc = 1, n_self = 4)
F1 <- as.numeric(EXP["HET"])
F2 <- as.numeric(EXP["ALT"])

# --- simulate BC1S4 ground truth on the experimental 118K v5 cM grid ----------
mc <- fread("data/teonam/markers_v5_gwas118k_cm.tsv")
setnames(mc, "pos_v5", "pos")
map <- data.frame(chr = as.integer(mc$chr), cm = as.numeric(mc$cm), bp = as.integer(mc$pos))
log_info("simulating %s n=%d on the %d-marker 118K grid (seed=%d) ...", DESIGN, N_RIL, nrow(map), SEED)
ts <- Sys.time()
truth_m <- as.data.table(simulate_nil(DESIGN, n = N_RIL, map = map, n_markers = nrow(map), seed = SEED))
truth <- as.data.table(to_segments(truth_m))
grid <- unique(truth_m[, .(chr = as.integer(chr), pos = as.integer(pos))])
setorder(grid, chr, pos)
setkey(truth_m, name, chr, pos)
log_info(
  "  truth: %d segments, %d markers, %d lines; %.1f breakpoints/line (%.1f s); priors f_1=%.3f f_2=%.3f",
  nrow(truth), nrow(grid), uniqueN(truth_m$name), breakpoint_count(truth) / uniqueN(truth_m$name),
  as.numeric(difftime(Sys.time(), ts, units = "secs")), F1, F2
)

# --- density-centered log grid for rrate (default rrate0 = 2*total_cM/(100*n_markers)) ---
total_cM <- sum(mc[, max(cm) - min(cm), by = chr]$V1)
rrate0 <- 2 * total_cM / (100 * nrow(mc))
RRATE_VALUES <- rrate0 * 2^(-4:3) # 8-point log2 grid bracketing the density default
log_info(
  "total_cM=%.0f, n_markers=%d -> rrate0=%.2e; rrate grid = %s",
  total_cM, nrow(mc), rrate0, paste(sprintf("%.2e", RRATE_VALUES), collapse = " ")
)
log_info(
  "=== nnil per-coverage rrate calibration: %d coverages x %d rrates + golden refine, %d threads ===",
  length(LAMBDAS), length(RRATE_VALUES), THREADS
)

# --- read draws for one coverage (matches the sweep's .draw_counts model) -----
make_counts <- function(lambda) {
  set.seed(SEED + 1000L * which(LAMBDAS == lambda))
  st <- truth_m$state
  if (is.infinite(lambda)) {
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

# --- per-coverage calibration: coarse sweep -> bracket -> golden refine --------
sweeps <- list()
best <- list()
t_loop <- Sys.time()
for (i in seq_along(LAMBDAS)) {
  lam <- LAMBDAS[i]
  covlab <- if (is.infinite(lam)) "Inf" else as.character(lam)
  d <- make_counts(lam)
  cov_frac <- mean(d$n_ref + d$n_alt > 0)
  tc <- Sys.time()
  log_info(
    "  [%d/%d] lambda=%-4s (%.0f%% cov): coarse sweep of %d rrates on %d threads ...",
    i, length(LAMBDAS), covlab, 100 * cov_frac, length(RRATE_VALUES), THREADS
  )
  sw <- sweep_calibrate(d, truth, grid,
    caller = "nnil", values = RRATE_VALUES,
    threads = THREADS, f_1 = F1, f_2 = F2, err = READ_PARS$error, min_reads = 1L
  )
  br <- bracket_from_sweep(sw, "donor_frag_FDR") # FDR is minimized
  log_info(
    "  [%d/%d] lambda=%-4s: coarse done (%.1f min); golden refine on rrate in [%.2e, %.2e] ...",
    i, length(LAMBDAS), covlab, as.numeric(difftime(Sys.time(), tc, units = "mins")),
    br[["lo"]], br[["hi"]]
  )
  ref <- golden_refine(d, truth, grid,
    caller = "nnil", lo = br[["lo"]], hi = br[["hi"]], objective = "donor_frag_FDR",
    threads = THREADS, f_1 = F1, f_2 = F2, err = READ_PARS$error, min_reads = 1L
  )
  el <- as.numeric(difftime(Sys.time(), tc, units = "mins"))
  sw[, coverage := covlab]
  sweeps[[covlab]] <- sw
  rrate_star <- ref$value
  fdr <- ref$score$donor_frag_FDR[1]
  rc <- ref$score$donor_marker_recall[1]
  dice <- ref$score$donor_frag_dice[1]
  best[[covlab]] <- data.table(
    coverage = covlab, cov_frac = round(cov_frac, 3),
    rrate = signif(rrate_star, 4), fdr = round(fdr, 4), recall = round(rc, 4),
    dice = round(dice, 4), mins = round(el, 1)
  )
  log_info(
    "  [%d/%d] lambda=%-4s (%.0f%% cov): rrate* = %.3e, FDR=%.3f recall=%.3f (%.1f min)",
    i, length(LAMBDAS), covlab, 100 * cov_frac, rrate_star, fdr, rc, el
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
  fwrite(B, file.path(OUT, "nnil_calib_bycov.csv"))
  fwrite(S, file.path(OUT, "nnil_calib_bycov_sweep.csv"))
  log_info("wrote nnil_calib_bycov.csv + _sweep.csv")
}
log_info("=== DONE: %d coverages in %.1f min total ===", length(LAMBDAS), as.numeric(difftime(Sys.time(), t_loop, units = "mins")))
print(B)
