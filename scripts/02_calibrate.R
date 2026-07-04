#!/usr/bin/env Rscript
# Heavy compute for analysis/02-simulation-calibration.qmd — kept OUT of the
# notebook so the render is fast (the note only reads the CSVs written here).
#
# Pipeline:
#   1. (optional) generate the BC2S2 truth + degraded skim counts   [R/simulate.R]
#      -> results/sim/skim.rds (compact bundle: grid once + count matrices) +
#         bc2s2_truth_segments.csv. Load with load_sim() / sim_counts().
#   2. two-stage calibrate each knob (nNIL rrate, RTIGER rigidity): log sweep
#      to bracket the optimum, then golden-ratio refine  [R/calibrate.R]
#   3. benchmark each caller at its refined F1-optimal knob against truth
#   4. write small summary tables + the fragment-size / genotype-frequency
#      tables the note plots.
#
# Outputs (results/sim/, gitignored — the rendered docs/ HTML carries the figures):
#   nnil_rrate_sweep.csv, rtiger_rigidity_sweep.csv   (stage-1 log-sweep curves)
#   nnil_rrate_refine.csv, rtiger_rigidity_refine.csv (stage-2 golden probes)
#   benchmark.csv, calib_params.csv, frag_sizes.csv, geno_fractions.csv
#
# Run:  Rscript scripts/02_calibrate.R            # uses existing results/sim
#       Rscript scripts/02_calibrate.R --generate # (re)generate the sim first

suppressMessages({
  library(nilHMM)
  library(data.table)
})
root <- here::here()
for (f in list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)) source(f)
SIM <- file.path(root, "results/sim")

# --- knobs (offline, so we can afford larger sets than a live render) --------
N_CAL_NNIL <- 300L # calibration subset for the nNIL log sweep
N_CAL_RTIG <- 60L # rtiger fit is joint EM -> smaller subset
N_BENCH <- 300L # held-out benchmark set
THREADS <- max(1L, parallel::detectCores() - 2L) # fan-out for caller_sweep
LOG_N <- 10L # coarse rrate log-sweep points before the golden refine
NNIL_VALUES <- log_grid(1e-6, 1e-1, LOG_N) # rrate: log-spaced (continuous)
RTIG_VALUES <- 2L^(1:9) # rigidity: powers of 2 (2..512); feasibility-filtered per cohort

if ("--generate" %in% commandArgs(TRUE)) {
  message("generating BC2S2 skim truth (n=1500, 50k markers)...")
  simulate_source("BC2S2", "skim", n = 1500L, seed = 1L)
}

sim <- fread(file.path(SIM, "bc2s2_truth_segments.csv")) # TRUTH segments
bundle <- load_sim(file.path(SIM, "skim.rds")) # compact counts bundle
grid <- as.data.table(bundle$grid)
ids <- sort(bundle$names)
# Materialize only the samples the sweeps/benchmark use (grid is shared, so a
# subset costs only its own columns); sub(skim, k) then filters within these.
NMAX <- max(N_CAL_NNIL, N_CAL_RTIG, N_BENCH)
skim <- as.data.table(sim_counts(bundle, ids[seq_len(min(NMAX, length(ids)))]))
message(sprintf(
  "loaded %d NILs (of %d), %d markers, %.0f%% missing",
  uniqueN(skim$name), length(ids), nrow(grid), 100 * mean(skim$n_ref + skim$n_alt == 0)
))

sub <- function(dt, k) dt[name %in% ids[seq_len(min(k, length(ids)))]]

# --- 1. two-stage calibration: log sweep -> golden-ratio refine --------------
# Stage 1 brackets the optimum on a coarse log grid (sweep_calibrate = one shared
# caller_sweep fit + truth scoring); stage 2 golden-section-refines within the
# bracket. Both maximize donor-fragment F1. [R/calibrate.R]
calibrate_param <- function(caller, cal_k, values, integer, ...) {
  sk <- sub(skim, cal_k)
  tr <- sub(sim, cal_k)
  # rtiger: keep only rigidities every chromosome can support (> 2*r covered markers)
  if (integer) values <- feasible_rigidity(sk, values)
  message(sprintf("  %s: log sweep (%d NILs, %d pts) ...", caller, uniqueN(sk$name), length(values)))
  sweep <- sweep_calibrate(sk, tr, grid,
    caller = caller, threads = THREADS, values = values, ...
  )
  br <- bracket_from_sweep(sweep, "donor_frag_F1")
  message(sprintf("  %s: golden refine in [%.3g, %.3g] ...", caller, br[["lo"]], br[["hi"]]))
  ref <- golden_refine(sk, tr, grid,
    caller = caller, lo = br[["lo"]], hi = br[["hi"]],
    objective = "donor_frag_F1", threads = THREADS, ...
  )
  list(sweep = sweep, refine = ref)
}
nnil_cal <- calibrate_param("nnil", N_CAL_NNIL, NNIL_VALUES, integer = FALSE, design = "BC2S2", err = 0.01)
rtig_cal <- calibrate_param("rtiger", N_CAL_RTIG, RTIG_VALUES, integer = TRUE, design = "BC2S2")
nnil_sweep <- nnil_cal$sweep
rtiger_sweep <- rtig_cal$sweep
rrate_star <- nnil_cal$refine$value
rig_star <- as.integer(rtig_cal$refine$value)
fwrite(nnil_sweep, file.path(SIM, "nnil_rrate_sweep.csv"))
fwrite(rtiger_sweep, file.path(SIM, "rtiger_rigidity_sweep.csv"))
fwrite(nnil_cal$refine$evals, file.path(SIM, "nnil_rrate_refine.csv"))
fwrite(rtig_cal$refine$evals, file.path(SIM, "rtiger_rigidity_refine.csv"))
message(sprintf("F1-optimal (refined): rrate=%.3g | rigidity=%d", rrate_star, rig_star))

# --- 2. benchmark at the calibrated knobs ------------------------------------
skim_b <- sub(skim, N_BENCH)
sim_b <- sub(sim, N_BENCH)
message("benchmark: calling ", N_BENCH, " NILs with each caller ...")
skim_nnil <- as.data.table(call_ancestry(skim_b, caller = "nnil", design = "BC2S2", rrate = rrate_star, err = 0.01))
skim_rtiger <- as.data.table(call_ancestry(skim_b, caller = "rtiger", design = "BC2S2", rigidity = rig_star))

score_caller <- function(called, tag, param) {
  mf <- marker_f1(called, sim_b, grid)
  ff <- donor_fragment_f1(called, sim_b)
  data.table(
    caller = tag, param = param,
    marker_macroF1 = round(mf$macro_f1, 3),
    donor_marker_F1 = round(mf$per_class[class == "donor(>0)", f1], 3),
    donor_marker_recall = round(mf$per_class[class == "donor(>0)", recall], 3),
    donor_frag_F1 = round(ff$f1, 3), donor_frag_FDR = round(ff$fdr, 3),
    ks_fragsize = round(fragment_size_ks(donor_block_sizes(called), donor_block_sizes(sim_b)), 3),
    breakpoints = breakpoint_count(called), breakpoints_truth = breakpoint_count(sim_b)
  )
}
benchmark <- rbindlist(list(
  score_caller(skim_nnil, "Skim-nNIL", rrate_star),
  score_caller(skim_rtiger, "Skim-RTIGER", rig_star)
))
fwrite(benchmark, file.path(SIM, "benchmark.csv"))
print(benchmark)

# --- 3. plot-ready tables (fragment sizes + genotype fractions) --------------
frag <- rbindlist(list(
  data.table(src = "truth", mb = donor_block_sizes(sim_b)),
  data.table(src = "Skim-nNIL", mb = donor_block_sizes(skim_nnil)),
  data.table(src = "Skim-RTIGER", mb = donor_block_sizes(skim_rtiger))
))
fwrite(frag, file.path(SIM, "frag_sizes.csv"))

geno <- rbindlist(list(
  cbind(src = "truth", genotype_fractions(sim_b)),
  cbind(src = "Skim-nNIL", genotype_fractions(skim_nnil)),
  cbind(src = "Skim-RTIGER", genotype_fractions(skim_rtiger))
))
fwrite(geno, file.path(SIM, "geno_fractions.csv"))

# --- 4. calibrated params + expectation for the note header ------------------
p0 <- single_locus_expectation(n_bc = 2L, n_self = 2L)
gf0 <- genotype_fractions(sim)
# data.frame, not data.table: `data.table(key = ...)` would treat `key` as the
# reserved key= argument (setkeyv error), not a column. The note reads $key/$value.
fwrite(data.frame(
  key = c(
    "rrate_star", "rigidity_star", "n_cal_nnil", "n_cal_rtiger", "n_bench",
    "REF_truth", "HET_truth", "ALT_truth", "dosage_truth",
    "REF_expect", "HET_expect", "ALT_expect", "dosage_expect"
  ),
  value = c(
    rrate_star, rig_star, N_CAL_NNIL, N_CAL_RTIG, N_BENCH,
    mean(gf0$REF), mean(gf0$HET), mean(gf0$ALT), mean(gf0$dosage),
    p0["REF"], p0["HET"], p0["ALT"], p0["dosage"]
  ), stringsAsFactors = FALSE
), file.path(SIM, "calib_params.csv"))

message("done. wrote sweep/benchmark/frag/geno/params CSVs to ", SIM)
