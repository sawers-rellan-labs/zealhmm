#!/usr/bin/env Rscript
# Heavy compute for analysis/02-simulation-calibration.qmd — kept OUT of the
# notebook so the render is fast (the note only reads the CSVs written here).
#
# Pipeline:
#   1. (optional) generate the BC2S2 truth + degraded skim counts   [R/simulate.R]
#      -> results/sim/skim.rds (compact bundle: grid once + count matrices) +
#         bc2s2_truth_segments.csv. Load with load_sim() / sim_counts().
#   2. sweep the duration knob for nNIL (rrate) and RTIGER (rigidity)
#   3. benchmark each caller at its F1-optimal knob against truth
#   4. write small summary tables + the fragment-size / genotype-frequency
#      tables the note plots.
#
# Outputs (results/sim/, gitignored — the rendered docs/ HTML carries the figures):
#   nnil_rrate_sweep.csv, rtiger_rigidity_sweep.csv
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
N_CAL_NNIL <- 300L # calibration subset for the (batched) nNIL sweep
N_CAL_RTIG <- 60L # rtiger EM is per-sample -> smaller sweep subset
N_BENCH <- 300L # held-out benchmark set
RRATE_GRID <- 10^seq(-6, -1.5, by = 0.5)
RIGIDITY_GRID <- c(20L, 50L, 100L, 200L, 400L) # min-run in markers (50k grid)

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

# --- 1. sweeps ---------------------------------------------------------------
message("nNIL rrate sweep on ", N_CAL_NNIL, " NILs ...")
nnil_sweep <- calibrate_sweep(sub(skim, N_CAL_NNIL), sub(sim, N_CAL_NNIL), grid,
  caller = "nnil", param = "rrate", values = RRATE_GRID, design = "BC2S2", err = 0.01
)
nnil_sweep[, truth_bp := breakpoint_count(sub(sim, N_CAL_NNIL))]
fwrite(nnil_sweep, file.path(SIM, "nnil_rrate_sweep.csv"))

message("RTIGER rigidity sweep on ", N_CAL_RTIG, " NILs ...")
rtiger_sweep <- calibrate_sweep(sub(skim, N_CAL_RTIG), sub(sim, N_CAL_RTIG), grid,
  caller = "rtiger", param = "rigidity", values = RIGIDITY_GRID, design = "BC2S2"
)
rtiger_sweep[, truth_bp := breakpoint_count(sub(sim, N_CAL_RTIG))]
fwrite(rtiger_sweep, file.path(SIM, "rtiger_rigidity_sweep.csv"))

rrate_star <- nnil_sweep[which.max(donor_frag_F1), value]
rig_star <- rtiger_sweep[which.max(donor_frag_F1), value]
message(sprintf("F1-optimal: rrate=%.2e | rigidity=%d", rrate_star, rig_star))

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
fwrite(data.table(
  key = c(
    "rrate_star", "rigidity_star", "n_cal_nnil", "n_cal_rtiger", "n_bench",
    "REF_truth", "HET_truth", "ALT_truth", "dosage_truth",
    "REF_expect", "HET_expect", "ALT_expect", "dosage_expect"
  ),
  value = c(
    rrate_star, rig_star, N_CAL_NNIL, N_CAL_RTIG, N_BENCH,
    mean(gf0$REF), mean(gf0$HET), mean(gf0$ALT), mean(gf0$dosage),
    p0["REF"], p0["HET"], p0["ALT"], p0["dosage"]
  )
), file.path(SIM, "calib_params.csv"))

message("done. wrote sweep/benchmark/frag/geno/params CSVs to ", SIM)
