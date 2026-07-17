#!/usr/bin/env Rscript
# Heavy compute for analysis/simulation-calibration.qmd — kept OUT of the
# notebook so the render is fast (the note only reads the CSVs written here).
#
# Pipeline:
#   1. (optional) generate the DESIGN truth + degraded skim counts   [R/simulate.R]
#      -> results/sim/skim.rds (compact bundle: grid once + count matrices) +
#         <design>_truth_segments.csv. Load with load_sim() / sim_counts().
#   2. two-stage calibrate each knob — nNIL rrate, RTIGER rigidity, LB-Impute
#      recombdist: log sweep to bracket the optimum, then golden-ratio refine
#      (maximize donor-fragment Dice)  [R/calibrate.R]
#   3. benchmark each caller at its refined Dice-optimal knob against truth
#   4. write small summary tables + the fragment-size / genotype-frequency tables
#
# DESIGN is parametrized (--design, default BC2S2). All design-derived quantities
# (start priors f_1/f_2, the hard-call breeding prior, expectations) are computed
# from it via breeding_prior() -- one function, any design, no lookup
# table. LB-Impute calibrates recombdist in cM (the sim lives on the consensus cM
# map), drp=TRUE (clean hom<->hom RIL breakpoints).
#
# Outputs (results/sim/, gitignored — the rendered docs/ HTML carries the figures):
#   nnil_rrate_sweep.csv, rtiger_rigidity_sweep.csv, lbimpute_recombdist_sweep.csv
#   *_refine.csv (golden probes), benchmark.csv, calib_params.csv,
#   frag_sizes.csv, geno_fractions.csv
#
# Run:  Rscript scripts/02_calibrate.R            # uses existing results/sim
#       Rscript scripts/02_calibrate.R --generate # (re)generate the sim first

suppressMessages({
  library(nilHMM)
  library(data.table)
})
root <- here::here()
for (f in list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)) source(f)
source(file.path(root, "scripts/logging.R")) # timestamped log_info + ETA idiom
SIM <- file.path(root, "results/sim")

# --- args: --design=X, --smoke (tiny fast run for ETA), --generate (rebuild sim) --
.args <- commandArgs(TRUE)
.getopt <- function(flag, default) {
  hit <- grep(sprintf("^%s=", flag), .args, value = TRUE)
  if (length(hit)) sub(sprintf("^%s=", flag), "", hit[1]) else default
}
# Single source of truth: the breeding design this whole script calibrates. Every
# design-derived quantity (simulated truth, start priors, the hard-call breeding
# prior, expectations) is COMPUTED from DESIGN via breeding_prior() -- no
# lookup tables, no divergence, works for any design. Parametrized so caller scripts
# can drive the runs: --design=BC2S3, --design=BC1S4, etc.
DESIGN <- .getopt("--design", "BC2S2")
SMOKE <- "--smoke" %in% .args
GENERATE <- "--generate" %in% .args
# smoke writes to a throwaway subdir so it never clobbers the real CSVs
OUT <- if (SMOKE) file.path(SIM, "smoke") else SIM
if (SMOKE) dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
.eta <- function(i, N, t0, tag) {
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  log_info(
    ">>> %s done (%d/%d) | elapsed %.1f min | avg %.1f min | ETA ~%.1f min remaining",
    tag, i, N, el, el / i, (el / i) * (N - i)
  )
}

# --- design + knobs (offline, so we can afford larger sets than a live render) --
design_prior <- breeding_prior(DESIGN) # design genotype-freq prior c(REF, HET, ALT)
F1 <- as.numeric(design_prior["HET"]) # HMM start prior (HET freq)
F2 <- as.numeric(design_prior["ALT"]) # HMM start prior (donor-hom freq)

# smoke = tiny sizes for a fast end-to-end check + per-caller timing basis
N_CAL_NNIL <- if (SMOKE) 20L else 300L # calibration subset for the geometric log sweep
N_CAL_RTIG <- if (SMOKE) 15L else 60L # rtiger fit is joint EM -> smaller subset
N_CAL_LBI <- if (SMOKE) 20L else 300L # lbimpute is exact per value (fast) -> large ok
N_BENCH <- if (SMOKE) 20L else 300L # held-out benchmark set
GEN_N <- if (SMOKE) 80L else 1500L # NILs to simulate under --generate
LOG_N <- if (SMOKE) 4L else 10L
THREADS <- max(1L, parallel::detectCores() - 2L)
NNIL_VALUES <- log_grid(1e-6, 1e-1, LOG_N) # rrate: log-spaced (continuous)
RTIG_VALUES <- 2L^(1:9) # rigidity: powers of 2; feasibility-filtered per cohort
LBI_VALUES <- log_grid(0.05, 50, if (SMOKE) 5L else 12L) # recombdist cM: [spacing, Haldane 50]

if (GENERATE) {
  log_info("generating %s skim truth (n=%d, 50k markers)%s ...", DESIGN, GEN_N, if (SMOKE) " [SMOKE]" else "")
  simulate_source(DESIGN, "skim", n = GEN_N, seed = 1L)
}

sim <- fread(file.path(SIM, sprintf("%s_truth_segments.csv", tolower(DESIGN)))) # TRUTH
bundle <- load_sim(file.path(SIM, "skim.rds")) # compact counts bundle
grid <- as.data.table(bundle$grid)
# cM per grid marker (sim lives on the consensus cM map) for the lbimpute unit="cm"
.map <- as.data.table(load_consensus_map())
.map[, chri := suppressWarnings(as.integer(gsub("\\D", "", as.character(chr))))]
to_cm <- bp_to_cm(.map[, .(chr = chri, bp, cm)]) # monotone bp -> cM Marey spline (per chr)
grid[, cm := to_cm(chr, pos)]
ids <- sort(bundle$names)
NMAX <- max(N_CAL_NNIL, N_CAL_RTIG, N_CAL_LBI, N_BENCH)
skim <- as.data.table(sim_counts(bundle, ids[seq_len(min(NMAX, length(ids)))]))
skim[grid, cm := i.cm, on = c("chr", "pos")] # lbimpute needs cM per marker
# --- hard-call for the categorical nNIL caller --------------------------------
# nNIL is categorical: it takes hard genotype calls, not read counts, so hard-call
# each marker to g in {0,1,2,3} (3 = missing). We pass the DESIGN's breeding prior
# explicitly -- call_gt's default is HWE, which is WRONG for NILs (no random mating:
# P(HET) is the scheme's expectation, not 2p(1-p)) and underperforms even a flat
# prior. Choosing the design prior is a conscious decision (not choosing = flat).
# Reuse design_prior -- the same genotype-freq prior that seeds F1/F2, computed once
# from breeding_prior(DESIGN); one source, no lookup table.
skim[, g := {
  gg <- call_gt(n_ref, n_alt, prior = design_prior, error = 0.01)
  gg[is.na(gg)] <- 3L # NA (zero depth) -> explicit missing state
  as.integer(gg)
}]
message(sprintf(
  "loaded %d NILs (of %d), %d markers, %.0f%% missing | design %s (f_1=%.3f f_2=%.3f)",
  uniqueN(skim$name), length(ids), nrow(grid), 100 * mean(skim$n_ref + skim$n_alt == 0),
  DESIGN, F1, F2
))

sub <- function(dt, k) dt[name %in% ids[seq_len(min(k, length(ids)))]]

# --- 1. two-stage calibration: log sweep -> golden-ratio refine --------------
calibrate_param <- function(caller, cal_k, values, integer, ...) {
  sk <- sub(skim, cal_k)
  tr <- sub(sim, cal_k)
  if (integer) values <- feasible_rigidity(sk, values) # rtiger only
  message(sprintf("  %s: log sweep (%d NILs, %d pts) ...", caller, uniqueN(sk$name), length(values)))
  sweep <- sweep_calibrate(sk, tr, grid,
    caller = caller, threads = THREADS, values = values, ...
  )
  br <- bracket_from_sweep(sweep, "donor_frag_dice")
  message(sprintf("  %s: golden refine in [%.3g, %.3g] ...", caller, br[["lo"]], br[["hi"]]))
  ref <- golden_refine(sk, tr, grid,
    caller = caller, lo = br[["lo"]], hi = br[["hi"]],
    objective = "donor_frag_dice", threads = THREADS, ...
  )
  list(sweep = sweep, refine = ref)
}
t0 <- Sys.time()
log_info("calibration start | design %s | 4 callers | smoke=%s | threads=%d", DESIGN, SMOKE, THREADS)
nnil_cal <- calibrate_param("nnil", N_CAL_NNIL, NNIL_VALUES, integer = FALSE, f_1 = F1, f_2 = F2, err = 0.01)
.eta(1L, 4L, t0, "nnil")
# bbNIL: same geometric rrate knob as nNIL, but count/BetaBinomial emission (the
# dispatcher-correct caller for skim read counts). Same grid/start-prior/err so
# the two differ only in emission (categorical hard-call vs BetaBinomial counts).
bbnil_cal <- calibrate_param("bbnil", N_CAL_NNIL, NNIL_VALUES, integer = FALSE, f_1 = F1, f_2 = F2, err = 0.01)
.eta(2L, 4L, t0, "bbnil")
rtig_cal <- calibrate_param("rtiger", N_CAL_RTIG, RTIG_VALUES, integer = TRUE)
.eta(3L, 4L, t0, "rtiger")
lbi_cal <- calibrate_param("lbimpute", N_CAL_LBI, LBI_VALUES,
  integer = FALSE, unit = "cm", err = 0.01, genotypeerr = 0.05, drp = TRUE
)
.eta(4L, 4L, t0, "lbimpute")
rrate_star <- nnil_cal$refine$value
rrate_bb_star <- bbnil_cal$refine$value
rig_star <- as.integer(rtig_cal$refine$value)
recomb_star <- lbi_cal$refine$value
fwrite(nnil_cal$sweep, file.path(OUT, "nnil_rrate_sweep.csv"))
fwrite(bbnil_cal$sweep, file.path(OUT, "bbnil_rrate_sweep.csv"))
fwrite(rtig_cal$sweep, file.path(OUT, "rtiger_rigidity_sweep.csv"))
fwrite(lbi_cal$sweep, file.path(OUT, "lbimpute_recombdist_sweep.csv"))
fwrite(nnil_cal$refine$evals, file.path(OUT, "nnil_rrate_refine.csv"))
fwrite(bbnil_cal$refine$evals, file.path(OUT, "bbnil_rrate_refine.csv"))
fwrite(rtig_cal$refine$evals, file.path(OUT, "rtiger_rigidity_refine.csv"))
fwrite(lbi_cal$refine$evals, file.path(OUT, "lbimpute_recombdist_refine.csv"))
message(sprintf(
  "Dice-optimal (refined): nnil rrate=%.3g | bbnil rrate=%.3g | rigidity=%d | recombdist=%.4g cM (drp=TRUE)",
  rrate_star, rrate_bb_star, rig_star, recomb_star
))

# --- 2. benchmark at the calibrated knobs ------------------------------------
skim_b <- sub(skim, N_BENCH)
sim_b <- sub(sim, N_BENCH)
message("benchmark: calling ", N_BENCH, " NILs with each caller ...")
skim_nnil <- as.data.table(call_ancestry(skim_b, caller = "nnil", f_1 = F1, f_2 = F2, rrate = rrate_star, err = 0.01))
skim_bbnil <- as.data.table(call_ancestry(skim_b, caller = "bbnil", f_1 = F1, f_2 = F2, rrate = rrate_bb_star, err = 0.01))
skim_rtiger <- as.data.table(call_ancestry(skim_b, caller = "rtiger", rigidity = rig_star))
skim_lbi <- as.data.table(call_ancestry(skim_b,
  caller = "lbimpute", unit = "cm",
  recombdist = recomb_star, err = 0.01, genotypeerr = 0.05, drp = TRUE
))

score_caller <- function(called, tag, param) {
  mf <- marker_dice(called, sim_b, grid)
  ff <- donor_fragment_dice(called, sim_b)
  data.table(
    caller = tag, param = param,
    marker_macro_dice = round(mf$macro_dice, 3),
    donor_marker_dice = round(mf$per_class[class == "donor(>0)", dice], 3),
    donor_marker_recall = round(mf$per_class[class == "donor(>0)", recall], 3),
    donor_frag_dice = round(ff$dice, 3), donor_frag_FDR = round(ff$fdr, 3),
    ks_fragsize = round(fragment_size_ks(donor_block_sizes(called), donor_block_sizes(sim_b)), 3),
    breakpoints = breakpoint_count(called), breakpoints_truth = breakpoint_count(sim_b)
  )
}
benchmark <- rbindlist(list(
  score_caller(skim_nnil, "Skim-nNIL", rrate_star),
  score_caller(skim_bbnil, "Skim-bbNIL", rrate_bb_star),
  score_caller(skim_rtiger, "Skim-RTIGER", rig_star),
  score_caller(skim_lbi, "Skim-LBimpute", recomb_star)
))
fwrite(benchmark, file.path(OUT, "benchmark.csv"))
print(benchmark)

# --- 3. plot-ready tables (fragment sizes + genotype fractions) --------------
frag <- rbindlist(list(
  data.table(src = "truth", mb = donor_block_sizes(sim_b)),
  data.table(src = "Skim-nNIL", mb = donor_block_sizes(skim_nnil)),
  data.table(src = "Skim-bbNIL", mb = donor_block_sizes(skim_bbnil)),
  data.table(src = "Skim-RTIGER", mb = donor_block_sizes(skim_rtiger)),
  data.table(src = "Skim-LBimpute", mb = donor_block_sizes(skim_lbi))
))
fwrite(frag, file.path(OUT, "frag_sizes.csv"))

geno <- rbindlist(list(
  cbind(src = "truth", genotype_fractions(sim)),
  cbind(src = "Skim-nNIL", genotype_fractions(skim_nnil)),
  cbind(src = "Skim-bbNIL", genotype_fractions(skim_bbnil)),
  cbind(src = "Skim-RTIGER", genotype_fractions(skim_rtiger)),
  cbind(src = "Skim-LBimpute", genotype_fractions(skim_lbi))
))
fwrite(geno, file.path(OUT, "geno_fractions.csv"))

# --- 4. calibrated params + expectation for the note header ------------------
p0 <- breeding_prior(DESIGN)
gf0 <- genotype_fractions(sim)
# data.frame, not data.table: `data.table(key = ...)` would treat `key` as the
# reserved key= argument. The note reads $key/$value.
cp_out <- data.frame(
  key = c(
    "design", "rrate_star", "rrate_bb_star", "rigidity_star", "recombdist_star", "lbimpute_drp",
    "n_cal_nnil", "n_cal_rtiger", "n_cal_lbimpute", "n_bench",
    "REF_truth", "HET_truth", "ALT_truth", "dosage_truth",
    "REF_expect", "HET_expect", "ALT_expect", "dosage_expect"
  ),
  value = c(
    DESIGN, rrate_star, rrate_bb_star, rig_star, recomb_star, "TRUE",
    N_CAL_NNIL, N_CAL_RTIG, N_CAL_LBI, N_BENCH,
    mean(gf0$REF), mean(gf0$HET), mean(gf0$ALT), mean(gf0$dosage),
    p0["REF"], p0["HET"], p0["ALT"], p0["ALT"] + p0["HET"] / 2
  ), stringsAsFactors = FALSE
)
# Two writes: the ACTIVE file consumers read, plus a DESIGN-STAMPED archive so a
# later calibration on a different design never silently clobbers this one. Every
# param file carries a `design` row, so provenance is always self-describing.
fwrite(cp_out, file.path(OUT, "calib_params.csv"))
fwrite(cp_out, file.path(OUT, sprintf("calib_params_%s.csv", tolower(DESIGN))))

log_info(
  "done. wrote calib_params.csv + calib_params_%s.csv (+ sweep/benchmark/frag/geno) to %s%s",
  tolower(DESIGN), OUT, if (SMOKE) " [SMOKE]" else ""
)
