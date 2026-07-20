#!/usr/bin/env Rscript
# Calibration foil, step 5: the SIM-calibrated operating point (NO chip).
#
# Simulate the nNIL breeding design on the EXACT nNIL marker set and genetic map,
# then calibrate the nnil caller's rrate against the DENSE simulated truth, and ask
# where the optimum lands relative to the chip-admissible band (step 4).
#
# Design      BC5S2 (n_bc=5, n_self=2) -- matches Holland's HMM priors f_1=1/128,
#             f_2~=0.0117 exactly (single_locus_expectation("BC5S2")); the physical
#             lines are labelled BC5F3/F4 but the caller models BC5S2 frequencies.
# Grid        the 63,904 lifted v5 nNIL markers with TeoNAM native cM (markers_v5.tsv);
#             identical density/positions to the real data, so rrate (per-adjacent-
#             marker) is directly comparable to the chip side.
# Degrade     observed hard genotypes g are drawn FORWARD through Holland's gt-
#             emission matrix at the FOUNDER-GENOTYPE non-informative rate
#             nir=0.594 (f0 measured on the chip NAM-founder genotypes,
#             agent/nnil_foil_estimate_nir.py; the union-ascertainment value, NOT
#             Holland's grid-selected 0.9; not a "biological" prescription) and
#             the real GBS missing rate. This reflects the true detectability of
#             donor segments; a count/depth degradation would misrepresent it
#             (donor allele != ALT: ~59% of donor markers are IBS with B73).
# Caller      nnil with Holland's emission (nir=0.9), SAME as the chip side, so the
#             only difference between the chip and sim calibration curves is the
#             truth source (sparse chip vs dense sim). This reproduces the real
#             situation: a chip-tuned caller applied to data whose true nir is 0.59.
#
#   Rscript scripts/nnil_foil/04_sim_calibrate.R
# Output (data/nnil_foil/):
#   sim_rrate_sweep.csv   rrate + donor_frag_dice/FDR, donor_marker_dice, macro_dice,
#                         n_breakpoints, ks_fragsize (vs dense sim truth)
#   sim_calib.json        rrate_sim*, single-locus check (expected vs observed), config
#   sim_truth_segments.csv, sim_frag_sizes.csv (for the step-6 figure)

suppressMessages({
  library(nilHMM)
  library(data.table)
  library(jsonlite)
})
root <- here::here()
for (f in list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)) source(f)
source(file.path(root, "scripts/logging.R"))
FOIL <- file.path(root, "data/nnil_foil")
EQUIV <- file.path(root, "data/nnil_equiv")

DESIGN <- "BC5S2"
NIR_FOUNDER <- 0.594 # f0 measured on the chip NAM-founder genotypes (sim's degradation nir)
N_SIM <- 1500L # NILs to simulate (user: keep 1500)
N_CAL <- 300L # subset used for the rrate sweep (plenty; matches 02_calibrate)
M_INT <- 10L # Stahl interference (maize-like), p=0
set.seed(1L)

# ---- marker grid (exact v5 nNIL markers + native cM) ------------------------
xw <- fread(file.path(FOIL, "markers_v5.tsv"))
markers <- data.table(chr = as.integer(xw$chr), bp = as.integer(xw$pos_v5), cm = xw$cm)[order(chr, bp)]
cmlen <- markers[, .(L = max(cm)), by = chr][order(chr)]
log_info("sim grid: %d markers, %d chr, total %.0f cM", nrow(markers), nrow(cmlen), sum(cmlen$L))

# ---- pedigree + Holland params ----------------------------------------------
pd <- parse_design(DESIGN)
bp <- .bcsft_pedigree(pd$n_bc, pd$n_self)
ped <- bp$ped
nid <- bp$nil_id
hp <- fromJSON(file.path(EQUIV, "params.json"))

# Holland's gt-emission matrix (3 true states x 4 obs {0,1,2,missing}), built at
# the founder-genotype nir for FORWARD generation of observed genotypes.
gt_emimat <- function(germ, gert, p, mr, nir) {
  matrix(c(
    (1 - germ) * (1 - mr), p * germ * (1 - mr), (1 - p) * germ * (1 - mr), mr,
    (((1 - nir) * 0.5 * gert) + nir * (1 - germ)) * (1 - mr),
    (((1 - nir) * (1 - gert)) + (nir * germ * p)) * (1 - mr),
    (((1 - nir) * 0.5 * gert) + nir * germ * (1 - p)) * (1 - mr), mr,
    ((1 - nir) * germ * (1 - p) + (nir * (1 - germ))) * (1 - mr), germ * p * (1 - mr),
    ((1 - nir) * (1 - germ) + (nir * germ * (1 - p))) * (1 - mr), mr
  ), nrow = 3, byrow = TRUE)
}
emimat <- gt_emimat(hp$germ, hp$gert, hp$p, hp$mr, NIR_FOUNDER)
draw_obs <- function(states) { # true dosage {0,1,2} -> observed g {0,1,2,3}
  g <- integer(length(states))
  for (s in 0:2) {
    ix <- which(states == s)
    if (length(ix)) g[ix] <- sample.int(4L, length(ix), replace = TRUE, prob = emimat[s + 1L, ]) - 1L
  }
  g
}

# ---- simulate N_SIM NILs: dense truth + forward-degraded observed g ---------
nms <- sprintf("sim%04d", seq_len(N_SIM))
truth_l <- vector("list", N_SIM)
data_l <- vector("list", N_SIM)
dose_tab <- integer(3) # pooled true-state counts, for the single-locus check
t0 <- Sys.time()
for (i in seq_len(N_SIM)) {
  dosage <- .simulate_dosage(ped, cmlen, markers, m = M_INT, p = 0, nil_id = nid)
  truth_l[[i]] <- .truth_segments(markers, dosage, nms[i])
  dose_tab <- dose_tab + tabulate(dosage + 1L, nbins = 3L)
  if (i <= N_CAL) {
    data_l[[i]] <- data.table(name = nms[i], chr = markers$chr, pos = markers$bp, g = draw_obs(dosage))
  }
  if (i %% 250 == 0 || i == N_SIM) {
    log_info("  simulated %d/%d (%.0fs)", i, N_SIM, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }
}
truth <- rbindlist(truth_l)
data <- rbindlist(data_l)
setorder(data, name, chr, pos)
fwrite(truth, file.path(FOIL, "sim_truth_segments.csv"))

# ---- single-locus check: does simcross reproduce BC5S2 expected freqs? ------
expd <- single_locus_expectation(DESIGN)
obsd <- dose_tab / sum(dose_tab)
log_info(
  "single-locus BC5S2 check | expected REF/HET/ALT = %.4f/%.4f/%.4f",
  expd["REF"], expd["HET"], expd["ALT"]
)
log_info(
  "                         | observed REF/HET/ALT = %.4f/%.4f/%.4f",
  obsd[1], obsd[2], obsd[3]
)
sl_ok <- max(abs(obsd - as.numeric(expd))) < 0.003
log_info("single-locus check %s (max abs diff %.4f)", if (sl_ok) "PASS" else "CHECK", max(abs(obsd - as.numeric(expd))))

# ---- calibrate rrate against the DENSE sim truth ----------------------------
grid_eval <- markers[, .(chr, pos = bp)] # full dense grid for marker_dice
cal_names <- nms[seq_len(N_CAL)]
tr_cal <- truth[name %in% cal_names]
values <- log_grid(1e-6, 1e-1, 24L)

score_one <- function(v) {
  called <- as.data.table(call_ancestry(
    data = data, caller = "nnil", rrate = v,
    germ = hp$germ, gert = hp$gert, p = hp$p, nir = hp$nir, mr = hp$mr,
    f_1 = hp$f_1, f_2 = hp$f_2
  ))
  mf <- marker_dice(called, tr_cal, grid_eval)
  ff <- donor_fragment_dice(called, tr_cal)
  dm <- mf$per_class[class == "donor(>0)"]
  data.table(
    rrate = v, donor_frag_dice = ff$dice, donor_frag_FDR = ff$fdr,
    donor_marker_dice = dm$dice, marker_macro_dice = mf$macro_dice,
    n_breakpoints = breakpoint_count(called),
    ks_fragsize = fragment_size_ks(donor_block_sizes(called), donor_block_sizes(tr_cal))
  )
}
t1 <- Sys.time()
log_info("sim-side rrate sweep: %d points on %d NILs (nnil, caller nir=0.9) ...", length(values), N_CAL)
sweep <- rbindlist(lapply(seq_along(values), function(i) {
  r <- score_one(values[i])
  log_info(
    "  rrate=%.3e | frag_dice=%.3f FDR=%.3f (%d/%d, %.0fs)",
    values[i], r$donor_frag_dice, r$donor_frag_FDR, i, length(values),
    as.numeric(difftime(Sys.time(), t1, units = "secs"))
  )
  r
}))
fwrite(sweep, file.path(FOIL, "sim_rrate_sweep.csv"))

rrate_sim <- sweep$rrate[which.max(sweep$donor_frag_dice)]
best_fd <- max(sweep$donor_frag_dice)
interior <- which.max(sweep$donor_frag_dice) %in% seq(2L, length(values) - 1L)
# donor_block_sizes() returns a numeric vector (Mb); wrap for fwrite
fwrite(data.table(block_mb = donor_block_sizes(truth)), file.path(FOIL, "sim_frag_sizes.csv"))

writeLines(toJSON(list(
  design = DESIGN, nir_founder = NIR_FOUNDER, nir_caller = hp$nir,
  n_sim = N_SIM, n_cal = N_CAL, n_markers = nrow(markers),
  rrate_sim_star = rrate_sim, frag_dice_best = best_fd, interior_optimum = interior,
  single_locus_expected = as.list(round(expd, 5)),
  single_locus_observed = list(REF = obsd[1], HET = obsd[2], ALT = obsd[3]),
  single_locus_pass = sl_ok
), auto_unbox = TRUE, digits = 8), file.path(FOIL, "sim_calib.json"))
log_info("rrate_sim* = %.4e (frag Dice %.3f, interior=%s)", rrate_sim, best_fd, interior)
