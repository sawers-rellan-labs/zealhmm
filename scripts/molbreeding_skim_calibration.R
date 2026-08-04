#!/usr/bin/env Rscript
# Holland-style nnil calibration for ZEAL, REAL data: MolBreeding target-seq (~110x)
# = TRUTH (Jim's chip analog), SNP50K skim (~0.4x) = TEST (Jim's GBS analog). Calibrate
# nnil against Jim's objective -- the per-cell/marker ancestry MISMATCH between the two
# platforms -- with FDR/recall/DSC/KS/breakpoints reported alongside (NOT as replacements).
#
# TRUTH  : nnil on the MolBreeding hard calls -> truth ancestry mosaic. The extracted 110x
#          genotypes g in {0,1,2} (extract_molbreeding_hardcalls.R) are the nnil INPUT (NOT
#          re-genotype-called); nnil smooths them to ancestry, mirroring Jim running the HMM on
#          the chip. Fixed config, r = molbreeding-panel map-derived ~1.67e-3.
# TEST   : skim counts -> ML genotype call (call_gt prior="flat" = argmax-GL; the BC2S3
#          design-prior MAP kept as a foil) -> nnil grid
#          (Holland emission grid crossed with a WIDE r grid, one point per order of
#          magnitude, to map r-sensitivity on the real data; SNP50K avg_r ~3.2e-4 is inside it).
# SCORE  : rasterize both onto the 9,157 informative sites; per (marker x line) cell,
#          pooled over the 14 calibration NIL pairs (calibration_pairing.csv, incl. the
#          PN4_SID330<->PN4_SID322 Jaccard pair). Primary = mismatch (1 - accuracy).
#
# Out: data/zeal/molbreeding_skim_calibration.csv  (one row per config; mismatch first)
#   Rscript scripts/molbreeding_skim_calibration.R
suppressMessages({
  library(devtools)
  library(data.table)
  library(parallel)
})
devtools::load_all(path.expand("~/repos/nilhmm"), quiet = TRUE) # caller_grid, call_gt, breeding_prior
root <- here::here()
for (f in list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)) source(f) # marker_dsc, donor_fragment_dsc, ...
source(file.path(root, "scripts/logging.R"))

N_CORES <- min(detectCores() - 2L, 8L)
DESIGN <- "BC2S3"
MB <- file.path(root, "data/zeal/molbreeding")
SK <- file.path(root, "data/zeal/skim")
R_MOLB <- 1.67e-3 # molbreeding-panel map-derived per-marker recombination (9,157 sites)
# Wide skim r grid (Fausto's intended points, one per order of magnitude) to map the
# r-sensitivity on the REAL data. The SNP50K map-derived avg_r ~3.2e-4 sits inside it;
# Holland calibrates r only in a narrow avg_r/2..avg_r*2 band, so this wide grid is the
# r-insensitivity diagnostic, not Holland's calibration band.
R_SKIM_SWEEP <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
TRUTH_CFG <- list(nir = 0.7, germ = 1e-3, gert = 1e-4, p = 0.9) # fixed truth caller config

# --- calibration pairs (truth molb sample -> test skim sample), 14 NIL pairs ---
pair <- fread(file.path(root, "data/zeal/correspondence/calibration_pairing.csv"))[in_calibration == TRUE]
log_info("calibration pairs: %d (truth molb <-> test skim; incl. Jaccard override)", nrow(pair))

# --- TRUTH: nnil on MolBreeding hard calls (fixed config), relabel to the test-sample name ---
hc <- fread(file.path(MB, "molbreeding_hardcalls_wsfilt.tsv")) # name, pedigree, marker, chr, pos, g
hc <- hc[name %in% pair$truth_sample]
hc[is.na(g), g := 3L]
grid_eval <- unique(hc[, .(chr, pos)])[order(chr, pos)] # 9,157 informative sites = evaluation grid
mr_molb <- 0 # 110x: ~0 structural missing (0.12% -> negligible)
truth_EG <- as.data.frame(c(TRUTH_CFG, list(mr = mr_molb)))
truth_seg <- as.data.table(caller_grid(hc[, .(name, chr, pos, g)],
  caller = "nnil",
  emission_grid = truth_EG, rrate = R_MOLB, design = DESIGN, threads = N_CORES
))
keep_seg <- c("name", "chr", "start_bp", "end_bp", "state")
truth <- truth_seg[, ..keep_seg]
truth[, name := pair$test_sample[match(name, pair$truth_sample)]] # relabel truth -> skim id
tr_sizes <- donor_block_sizes(truth)
log_info(
  "TRUTH built: nnil on molbreeding hard calls | cfg nir=%.2f germ=%.0e gert=%.0e p=%.2f r=%.2e | %d samples",
  TRUTH_CFG$nir, TRUTH_CFG$germ, TRUTH_CFG$gert, TRUTH_CFG$p, R_MOLB, uniqueN(truth$name)
)

# --- TEST: skim counts (prior-independent load); hard call compared under two priors ---
skim0 <- rbindlist(lapply(pair$test_sample, function(s) {
  cf <- fread(file.path(SK, "counts", paste0(s, ".tsv")),
    header = FALSE,
    col.names = c("contig", "pos", "rb", "rc", "ab", "ac")
  )
  data.table(
    name = s, chr = as.integer(sub("chr", "", cf$contig)), pos = cf$pos,
    n_ref = as.integer(cf$rc), n_alt = as.integer(cf$ac)
  )
}))
mr_skim <- skim0[, mean(n_ref + n_alt == 0)]
design_prior <- breeding_prior(DESIGN)
# Skim genotype caller. ML = call_gt(prior = "flat") = uniform-prior argmax-GL = the
# maximum-likelihood genotype call, the standard data-only caller. ML IS THE CALIBRATION
# BASIS, for three documented reasons (the head-to-head vs `breeding` below substantiates 2-3):
#   (1) FAITHFUL to Holland's pipeline -- Jim called GBS genotypes with a standard caller,
#       NOT a breeding-design prior;
#   (2) it DOMINATES the BC2S3 breeding-design prior on every metric (mismatch, donor
#       recall, DSC, FDR, KS) and recovers nir = 0.7, matching the MolBreeding truth caller
#       and the authentic donor non-informative rate (~0.68);
#   (3) the breeding prior's het-suppression drives het_recall to 0 -- it bakes the expected
#       BC2S3 mosaic into the TEST before calibration, which is circular.
# `breeding` = the BC2S3 breeding-design prior, kept only as the het-suppressing foil that
# documents reasons (2)-(3).
PRIORS <- list(ml = "flat", breeding = design_prior)
log_info(
  "TEST: %d skim libraries, %d markers | realized skim missing mr=%.3f | comparing priors: %s",
  uniqueN(skim0$name), uniqueN(skim0[, .(chr, pos)]), mr_skim, paste(names(PRIORS), collapse = ", ")
)

# --- Holland emission grid x wide r sweep ---
EG <- as.data.table(expand.grid(
  nir = c(0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9),
  germ = c(1e-4, 1e-3, 1e-2), gert = c(1e-4, 1e-3, 1e-2),
  p = c(0.1, 0.25, 0.5, 0.75, 0.9)
))
EG[, mr := mr_skim]
truth_bp <- breakpoint_count(truth)
BATCH <- 21L
batches <- split(seq_len(nrow(EG)), (seq_len(nrow(EG)) - 1L) %/% BATCH)

# run the full grid for one skim hard-call prior; PRIMARY metric = per-cell mismatch
run_prior <- function(prior_vec, prior_name) {
  sk <- copy(skim0)
  sk[, g := {
    gg <- call_gt(n_ref, n_alt, prior = prior_vec, error = 0.01)
    gg[is.na(gg)] <- 3L
    as.integer(gg)
  }]
  setorder(sk, name, chr, pos)
  gdat <- sk[, .(name, chr, pos, g)]
  cov <- gdat[g != 3L]
  log_info(
    "[%s prior] skim hard-call composition (of covered): REF %.3f HET %.3f ALT %.3f",
    prior_name, mean(cov$g == 0), mean(cov$g == 1), mean(cov$g == 2)
  )
  score_batch <- function(eg_batch) {
    seg <- as.data.table(caller_grid(gdat,
      caller = "nnil", emission_grid = as.data.frame(eg_batch),
      rrate = R_SKIM_SWEEP, design = DESIGN, threads = N_CORES
    ))
    seg[, cfg := paste(nir, germ, gert, p, rrate, sep = "_")]
    cfgs <- unique(seg[, .(nir, germ, gert, p, r = rrate, cfg)])
    rbindlist(lapply(seq_len(nrow(cfgs)), function(i) {
      called <- seg[cfg == cfgs$cfg[i]]
      md <- marker_dsc(called, truth, grid_eval)
      pc <- md$per_class
      ff <- donor_fragment_dsc(called, truth)
      data.table(
        prior = prior_name,
        nir = cfgs$nir[i], germ = cfgs$germ[i], gert = cfgs$gert[i], p = cfgs$p[i], r = cfgs$r[i],
        mismatch = 1 - md$accuracy, # PRIMARY: Holland skim-vs-molbreeding per-cell mismatch
        donor_marker_recall = pc[class == "donor(>0)"]$recall,
        het_recall = pc[class == "HET"]$recall,
        alt_recall = pc[class == "ALT"]$recall,
        donor_frag_dsc = ff$dsc, donor_frag_FDR = ff$fdr,
        ks_fragsize = fragment_size_ks(donor_block_sizes(called), tr_sizes),
        breakpoint_ratio = breakpoint_count(called) / truth_bp
      )
    }))
  }
  t0 <- Sys.time()
  res <- vector("list", length(batches))
  for (b in seq_along(batches)) {
    res[[b]] <- score_batch(EG[batches[[b]]])
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    log_info(
      "[%s] batch %d/%d | elapsed %.1f min | ETA ~%.1f min", prior_name, b, length(batches), el,
      (el / b) * (length(batches) - b)
    )
  }
  rbindlist(res)
}

score <- rbindlist(lapply(names(PRIORS), function(nm) run_prior(PRIORS[[nm]], nm)))
setorder(score, prior, mismatch)
fwrite(score, file.path(root, "data/zeal/molbreeding_skim_calibration.csv"))
log_info(
  "wrote data/zeal/molbreeding_skim_calibration.csv (%d rows = %d configs x %d priors)",
  nrow(score), nrow(EG) * length(R_SKIM_SWEEP), length(PRIORS)
)

# --- comparison: design vs flat skim hard-call prior ---
cat("\n=== BEST by mismatch, per prior ===\n")
print(score[, .SD[which.min(mismatch)], by = prior][, .(prior,
  mismatch = round(mismatch, 4),
  nir, p, r, het_recall = round(het_recall, 3), alt_recall = round(alt_recall, 3),
  donor_marker_recall = round(donor_marker_recall, 3), donor_frag_FDR = round(donor_frag_FDR, 3)
)])
cat("\n=== marginal mismatch & het_recall vs nir, by prior ===\n")
print(dcast(score[, .(mismatch = mean(mismatch), het_recall = mean(het_recall)), by = .(prior, nir)],
  nir ~ prior,
  value.var = c("mismatch", "het_recall")
))
