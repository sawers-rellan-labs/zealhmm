#!/usr/bin/env Rscript
# =============================================================================
# Per-coverage LB-Impute recombdist calibration on SIMULATED ground truth (118K grid).
#
# The LB-Impute analog of scripts/teonam_{rtiger,nnil}_calib_bycov.R. The tuned knob
# is `recombdist` (cM): the linkage-decay DISTANCE over which the distance-dependent
# transition relaxes toward the stationary switch probability. Small recombdist ->
# the transition relaxes quickly -> shorter/more segments; large recombdist ->
# stiffer -> longer segments. Unlike RTIGER's index-based rigidity and nNIL's
# coordinate-free rrate, LB-Impute is the MAP-AWARE caller (unit = "cm" consumes the
# genetic map's cM gaps), so its parameter lives on the cM scale and its grid is
# centred on the ~50 cM default. recombdist is CONTINUOUS and touches only the
# transition (never the emission), so every swept value is EXACT (identical to a cold
# call_ancestry at that recombdist) and we keep the golden-ratio refine.
#
# Ground truth is a SIMULATION (simulate_nil, known crossovers) on the SAME 118K grid
# and read model (.draw_counts, error=0.01) as the sweep; objective = MIN donor-fragment
# FDR (err-long). LB-Impute keeps zero-read markers (flat emission) so the cM-distance
# transition sees true marker spacing -> min_reads is a NO-OP (pass 0L; all markers kept).
#
# Output (ephemeral sim, only small CSVs written):
#   results/sim/teonam/lbimpute_calib_bycov.csv        (recombdist*(lambda) + FDR/recall)
#   results/sim/teonam/lbimpute_calib_bycov_sweep.csv  (full coarse FDR curves)
# Run: Rscript scripts/teonam_lbimpute_calib_bycov.R [--smoke]   (smoke = lambda=1 only)
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
# LB-Impute emission/transition parameters (fixed across the recombdist grid).
GENOTYPEERR <- 0.05 # emission floor/ceiling (errg)
DRP <- TRUE # double-recombination penalty
# caller_sweep segments the lbimpute states per-column (n_markers x V matrix), so the
# master holds ~24M rows, not the n_values-stacked matrix that OOMed nNIL pre-fix.
# Same RTIGER thread convention; env-overridable if memory needs capping.
THREADS <- as.integer(Sys.getenv("LBIMPUTE_THREADS", as.character(max(1L, parallel::detectCores() - 2L))))

# LB-Impute start-distribution seed for BC1S4 (the Mendelian breeding_prior).
# For lbimpute f_1/f_2 only SEED the start dist (no state-frequency prior in the
# transition) -- same source as the nNIL/rtiger calibrations.
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
# NB: simulate_nil() resamples an evenly-spaced marker grid (build_marker_grid, hyman
# bp->cM spline) and RETURNS a `cm` column per marker -- it does NOT reuse mc$pos, so a
# (chr,pos) join back to mc misses (NA) and a cm join would be many-to-many on tied cM.
# Carry the sim's own cm (the native map the truth was drawn on) row-aligned into the
# per-coverage counts below. Same truth_m (seed 12345) as the rtiger/nnil calibrations.
log_info(
  "  truth: %d segments, %d markers, %d lines; %.1f breakpoints/line (%.1f s); start-seed f_1=%.3f f_2=%.3f",
  nrow(truth), nrow(grid), uniqueN(truth_m$name), breakpoint_count(truth) / uniqueN(truth_m$name),
  as.numeric(difftime(Sys.time(), ts, units = "secs")), F1, F2
)

# --- cM-scale log grid for recombdist, centred on the ~50 cM default ----------
# Geometric centre sqrt(8*256) ~ 45 cM (near the ~50 default); span brackets both a
# fast-relaxing (8 cM) and a stiff (256 cM) transition. larger recombdist = stiffer.
RECOMBDIST_VALUES <- log_grid(8, 256, 8L)
log_info(
  "recombdist grid (cM) = %s",
  paste(sprintf("%.1f", RECOMBDIST_VALUES), collapse = " ")
)
log_info(
  "=== LB-Impute per-coverage recombdist calibration: %d coverages x %d recombdists + golden refine, %d threads ===",
  length(LAMBDAS), length(RECOMBDIST_VALUES), THREADS
)

# --- read draws for one coverage (matches the sweep's .draw_counts model) ------
# Keeps ALL markers (incl zero-read) + attaches cm: lbimpute uses the flat emission
# on zero-read markers so the cM transition sees true spacing (min_reads no-op).
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
    n_ref = n_ref, n_alt = n_alt, cm = as.numeric(truth_m$cm) # sim's native-map cM, row-aligned
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
    "  [%d/%d] lambda=%-4s (%.0f%% cov): coarse sweep of %d recombdists on %d threads ...",
    i, length(LAMBDAS), covlab, 100 * cov_frac, length(RECOMBDIST_VALUES), THREADS
  )
  sw <- sweep_calibrate(d, truth, grid,
    caller = "lbimpute", values = RECOMBDIST_VALUES, threads = THREADS,
    unit = "cm", genotypeerr = GENOTYPEERR, drp = DRP, err = READ_PARS$error,
    f_1 = F1, f_2 = F2, min_reads = 0L
  )
  br <- bracket_from_sweep(sw, "donor_frag_FDR") # FDR is minimized
  log_info(
    "  [%d/%d] lambda=%-4s: coarse done (%.1f min); golden refine on recombdist in [%.1f, %.1f] cM ...",
    i, length(LAMBDAS), covlab, as.numeric(difftime(Sys.time(), tc, units = "mins")),
    br[["lo"]], br[["hi"]]
  )
  ref <- golden_refine(d, truth, grid,
    caller = "lbimpute", lo = br[["lo"]], hi = br[["hi"]], objective = "donor_frag_FDR",
    threads = THREADS, unit = "cm", genotypeerr = GENOTYPEERR, drp = DRP,
    err = READ_PARS$error, f_1 = F1, f_2 = F2, min_reads = 0L
  )
  el <- as.numeric(difftime(Sys.time(), tc, units = "mins"))
  sw[, coverage := covlab]
  sweeps[[covlab]] <- sw
  rd_star <- ref$value
  fdr <- ref$score$donor_frag_FDR[1]
  rc <- ref$score$donor_marker_recall[1]
  dice <- ref$score$donor_frag_dice[1]
  best[[covlab]] <- data.table(
    coverage = covlab, cov_frac = round(cov_frac, 3),
    recombdist = signif(rd_star, 4), fdr = round(fdr, 4), recall = round(rc, 4),
    dice = round(dice, 4), mins = round(el, 1)
  )
  log_info(
    "  [%d/%d] lambda=%-4s (%.0f%% cov): recombdist* = %.2f cM, FDR=%.3f recall=%.3f (%.1f min)",
    i, length(LAMBDAS), covlab, 100 * cov_frac, rd_star, fdr, rc, el
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
  fwrite(B, file.path(OUT, "lbimpute_calib_bycov.csv"))
  fwrite(S, file.path(OUT, "lbimpute_calib_bycov_sweep.csv"))
  log_info("wrote lbimpute_calib_bycov.csv + _sweep.csv")
}
log_info("=== DONE: %d coverages in %.1f min total ===", length(LAMBDAS), as.numeric(difftime(Sys.time(), t_loop, units = "mins")))
print(B)
