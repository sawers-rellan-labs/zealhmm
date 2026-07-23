#!/usr/bin/env Rscript
# ZEAL molbreeding-skim nnil calibration figure (analog of fig_nnil_sim_calibration.R).
#
# Two calibration truths:
#   SIM  = the BC2S3 single-individual baseline (scripts/simulate_zeal_nil.R); truth =
#          the latent ancestry mosaic (dose 0/1/2). nnil runs on ML-hard-called sim counts.
#   MolB = nnil calls on the MolBreeding hard-call genotypes (the real independent truth).
# Test = the 14 real skim NILs (calibration_pairing.csv), ML-hard-called (flat prior).
# ML (flat) hard-call is applied identically to sim counts and skim counts.
#
# ZEAL has TWO calibration points because the criteria disagree:
#   mismatch-optimum  nir* = argmin marker mismatch (Holland's criterion)
#   fragment-optimum  nir* = argmax donor-fragment DSC
# Both operating points are carried through panels B/C/D.
#
# Panels (2x2):
#   A  marker mismatch vs nir, nnil scored against the SIM truth and the MolB truth;
#      open circles = each truth's mismatch-opt;
#      dashed vertical = taxon nir ~0.69 (mean REF+0.5*HET over all taxa).
#   B  SIM: introgression-size ECDF, nnil at the two operating points vs the sim latent truth.
#   C  Real skim: introgression-size ECDF, nnil at the two operating points vs the MolB truth.
#   D  QQ of introgression size, mismatch-cal vs fragment-cal, on the real skim.
#
# Compute is CACHED; NNIL_ZEALCAL_RECOMPUTE=1 forces a rebuild.
#   Rscript scripts/fig_zeal_molbreeding_skim_calibration.R

suppressMessages({
  library(devtools)
  library(data.table)
  library(here)
  library(ggplot2)
  library(patchwork)
})
devtools::load_all(path.expand("~/repos/nilhmm"), quiet = TRUE) # caller_grid, call_gt (feat branch, not installed nilHMM)
source(here::here("R/metrics.R"))
source(here::here("scripts/logging.R"))

OUT <- here::here("results/sim/zeal_nil")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
NIR_GRID <- c(0.01, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
FIX <- list(germ = 1e-4, gert = 1e-2, p = 0.1) # fixed emission background (calibration mismatch-opt non-nir); sweep nir
R_MOLB <- 1.67e-3
TRUTH_CFG <- data.frame(nir = 0.7, germ = 1e-3, gert = 1e-4, p = 0.9, mr = 0) # MolB truth caller
KEEP <- c("name", "chr", "start_bp", "end_bp", "state")
BASE <- 20
CACHE <- file.path(OUT, "fig_zeal_molbreeding_skim_calibration_cache.rds")

if (file.exists(CACHE) && Sys.getenv("NNIL_ZEALCAL_RECOMPUTE") == "") {
  log_info("[fig] reusing cached compute (%s)", basename(CACHE))
  D <- readRDS(CACHE)
} else {
  ml_g <- function(nr, na) {
    g <- nilHMM::call_gt(as.integer(nr), as.integer(na), prior = "flat", error = 0.01, return = "call")
    g[is.na(g)] <- 3L
    as.integer(g)
  }
  seg_at <- function(seg, v) seg[abs(nir - v) < 1e-9]

  # ---- SIM leg: ML-called counts -> nnil sweep; truth = latent ancestry ----
  # Subsample lines for the calibration curve (a pooled average; a few hundred lines
  # give the same mismatch/DSC shape). NNIL_ZEALCAL_NSIM overrides (default 500).
  sim <- readRDS(file.path(OUT, "zeal_nil_bc2s3_full.rds"))
  N_SIM <- as.integer(Sys.getenv("NNIL_ZEALCAL_NSIM", "500"))
  set.seed(1)
  sub <- sort(sample(seq_along(sim$names), min(N_SIM, length(sim$names))))
  sub_names <- sim$names[sub]
  M <- nrow(sim$grid)
  sim_long <- data.table(
    name = rep(sub_names, each = M),
    chr = rep(as.integer(sim$grid$chr), length(sub)),
    pos = rep(as.integer(sim$grid$pos), length(sub)),
    g = ml_g(sim$n_ref[, sub], sim$n_alt[, sub])
  )
  sim_truth <- as.data.table(sim$truth)[name %in% sub_names]
  sim_grid <- data.table(chr = as.integer(sim$grid$chr), pos = as.integer(sim$grid$pos))
  sim_traster <- rasterize_named(sim_truth, sim_grid) # rasterize the (constant) sim truth ONCE
  R_SKIM <- sim$map_r # SNP50K map-derived r (same cM map as the skim)
  EG <- data.frame(
    nir = NIR_GRID, germ = FIX$germ, gert = FIX$gert, p = FIX$p,
    mr = mean((sim$n_ref[, sub] + sim$n_alt[, sub]) == 0)
  )
  log_info("[fig] SIM: %d/%d lines x %d markers; r=%.3e; nir grid %d", length(sub), length(sim$names), M, R_SKIM, length(NIR_GRID))
  threads <- min(parallel::detectCores() - 2L, 8L)
  t0 <- Sys.time()
  sim_seg <- rbindlist(lapply(seq_along(NIR_GRID), function(i) {
    s <- as.data.table(caller_grid(sim_long,
      caller = "nnil", emission_grid = EG[i, , drop = FALSE],
      rrate = R_SKIM, design = "BC2S3", threads = threads
    ))
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    log_info(
      "[fig] sim decode nir=%.2f (%d/%d) | elapsed %.0fs | ETA ~%.0fs",
      NIR_GRID[i], i, length(NIR_GRID), el, el / i * (length(NIR_GRID) - i)
    )
    s
  }))
  sim_mm <- sapply(NIR_GRID, function(v) 1 - marker_dice(seg_at(sim_seg, v)[, ..KEEP], sim_truth, sim_grid, truth_raster = sim_traster)$accuracy)

  # ---- MolB truth: nnil on molbreeding hard calls; skim test: ML-called counts ----
  pair <- fread(here::here("data/zeal/correspondence/calibration_pairing.csv"))[in_calibration == TRUE]
  hc <- fread(here::here("data/zeal/molbreeding/molbreeding_hardcalls_wsfilt.tsv"))[name %in% pair$truth_sample]
  hc[is.na(g), g := 3L]
  molb_grid <- unique(hc[, .(chr, pos)])[order(chr, pos)]
  molb_truth <- as.data.table(caller_grid(hc[, .(name, chr, pos, g)],
    caller = "nnil",
    emission_grid = TRUTH_CFG, rrate = R_MOLB, design = "BC2S3", threads = 4L
  ))[, ..KEEP]
  molb_truth[, name := pair$test_sample[match(name, pair$truth_sample)]]

  skim <- rbindlist(lapply(pair$test_sample, function(s) {
    cf <- fread(here::here("data/zeal/skim/counts", paste0(s, ".tsv")),
      header = FALSE,
      col.names = c("contig", "pos", "rb", "rc", "ab", "ac")
    )
    data.table(
      name = s, chr = as.integer(sub("chr", "", cf$contig)), pos = cf$pos,
      g = ml_g(cf$rc, cf$ac)
    )
  }))
  setorder(skim, name, chr, pos)
  EGk <- data.frame(
    nir = NIR_GRID, germ = FIX$germ, gert = FIX$gert, p = FIX$p,
    mr = mean(skim$g == 3L)
  )
  skim_seg <- as.data.table(caller_grid(skim[, .(name, chr, pos, g)],
    caller = "nnil",
    emission_grid = EGk, rrate = R_SKIM, design = "BC2S3", threads = 4L
  ))
  molb_traster <- rasterize_named(molb_truth, molb_grid) # rasterize the MolB truth ONCE
  molb_tblocks <- .donor_blocks(molb_truth) # merge MolB truth blocks ONCE
  skim_mm <- sapply(NIR_GRID, function(v) 1 - marker_dice(seg_at(skim_seg, v)[, ..KEEP], molb_truth, molb_grid, truth_raster = molb_traster)$accuracy)
  skim_dsc <- sapply(NIR_GRID, function(v) donor_fragment_dice(seg_at(skim_seg, v)[, ..KEEP], molb_truth, truth_blocks = molb_tblocks)$dice)

  # Two calibrations (nnil-figure analog): calibrate on the SIM and on MolBreeding, then
  # apply BOTH to the real skim (panels C/D). nir_frag (MolBreeding DSC-opt) is marked in A only.
  nir_sim <- NIR_GRID[which.min(sim_mm)]
  nir_molb <- NIR_GRID[which.min(skim_mm)]
  nir_frag <- NIR_GRID[which.max(skim_dsc)]
  D <- list(
    sweep = rbind(
      data.table(truth = "sim", nir = NIR_GRID, mismatch = sim_mm),
      data.table(truth = "molb", nir = NIR_GRID, mismatch = skim_mm)
    ),
    nir_sim = nir_sim, nir_molb = nir_molb, nir_frag = nir_frag,
    sim_truth_sz = donor_block_sizes(sim_truth),
    sim_sz_sim = donor_block_sizes(seg_at(sim_seg, nir_sim)[, ..KEEP]), # sim @ sim-cal
    sim_sz_molb = donor_block_sizes(seg_at(sim_seg, nir_molb)[, ..KEEP]), # sim @ MolB-cal
    molb_truth_sz = donor_block_sizes(molb_truth),
    skim_sz_sim = donor_block_sizes(seg_at(skim_seg, nir_sim)[, ..KEEP]), # real skim @ sim-cal
    skim_sz_molb = donor_block_sizes(seg_at(skim_seg, nir_molb)[, ..KEEP]), # real skim @ MolB-cal
    nir_donor = 0.676 # authentic donor non-informative rate (f_REF + 0.5 f_HET), ZEAL cohort
  )
  saveRDS(D, CACHE)
  fwrite(D$sweep, file.path(OUT, "fig_zeal_molbreeding_skim_calibration_sweep.csv"))
}
list2env(D, environment())
log_info("[fig] optima: sim mismatch nir=%.2f | MolBreeding mismatch nir=%.2f | MolBreeding DSC nir=%.2f", nir_sim, nir_molb, nir_frag)

# taxon nir = mean of (f_REF + 0.5*f_HET) over ALL taxa, from the reference-panel per-taxon
# genotype frequencies (grand-pooled over every panel individual of the 5 sim taxa). This is the
# non-informative-marker rate of a taxon-level donor, the ZEAL analog of the nNIL per-donor mean
# REF fraction; the mismatch/DSC optima are read against this single reference.
.gp <- readRDS(here::here("data/zeal/reference_panel/ref_panel_gt.rds"))
.rid <- fread(here::here("data/zeal/reference_panel/reference_ids_by_taxon.csv"))
.samps <- intersect(.rid[maizegdb_prefix %in% c("Zx", "Zv", "Zd", "Zl", "Zh")]$sample, colnames(.gp$mat))
.sub <- .gp$mat[, .samps, drop = FALSE]
.ncov <- rowSums(.sub != 3L)
nir_taxon <- mean(((rowSums(.sub == 0L) + 0.5 * rowSums(.sub == 1L)) / .ncov)[.ncov > 0])
log_info("[fig] taxon nir (REF+0.5*HET, all-taxa reference-panel genotype freqs) = %.3f", nir_taxon)

# ================================ plot ======================================
# colour is CONSISTENT across all panels: blue = sim (calibration source / sim-calibrated
# caller), orange = MolBreeding (source / MolBreeding-calibrated caller), black dotted = truth.
col_sim <- "#0072B2"
col_molb <- "#D55E00"
L_sim <- sprintf("sim calibrated \nnnil (nir=%.2f)", nir_sim)
L_molb <- sprintf("molb calibrated \nnnil (nir=%.2f)", nir_molb)
L_st <- "sim latent ancestry\nBC2S3"
L_mt <- "molb calls"
pal <- c(setNames(c(col_sim, col_molb), c(L_sim, L_molb)), setNames(c("black", "black"), c(L_st, L_mt)))
lty <- c(setNames(c("solid", "solid"), c(L_sim, L_molb)), setNames(c("dotted", "dotted"), c(L_st, L_mt)))
MB_BREAKS <- c(0.1, 1, 10, 100)
allsz <- c(sim_truth_sz, sim_sz_sim, sim_sz_molb, molb_truth_sz, skim_sz_sim, skim_sz_molb)
xlim_mb <- c(0.1, max(allsz[is.finite(allsz) & allsz > 0]))

# A: sim (blue) + MolBreeding (orange) mismatch curves; circles = mismatch-opt; triangle = MolB DSC-opt
opt <- sweep[, .SD[which.min(mismatch)], by = truth]
p_A <- ggplot(sweep, aes(nir, mismatch, colour = truth)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.4) +
  geom_point(data = opt, size = 3.6, shape = 21, fill = "white", stroke = 1.1) +
  geom_vline(xintercept = nir_taxon, linetype = "dashed", colour = "black", linewidth = 0.5, alpha = 0.45) +
  annotate("text",
    x = nir_taxon, y = 0.30, label = "taxon",
    angle = 90, hjust = 0.5, vjust = 1.4, size = BASE * 0.24, colour = "black", alpha = 0.45
  ) +
  scale_colour_manual(
    values = c(sim = col_sim, molb = col_molb),
    labels = c(sim = "sim skim vs sim ancestry", molb = "real skim vs molb ancestry"), name = "calibration"
  ) +
  labs(x = expression("non-informative rate " * italic(nir)), y = "marker mismatch rate", title = "nnil nir calibration") +
  theme_bw(base_size = BASE) +
  theme(
    aspect.ratio = 1, legend.position = c(0.5, 0.99), legend.justification = c(0.5, 1),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )

mk_ecdf <- function(d, ttl, lvls) {
  d[, series := factor(series, levels = lvls)]
  ggplot(d, aes(size_mb, colour = series, linetype = series)) +
    stat_ecdf(linewidth = 1) +
    scale_x_log10(breaks = MB_BREAKS, limits = xlim_mb, oob = scales::oob_keep) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_linetype_manual(values = lty, name = NULL) +
    labs(x = "introgression size (Mb)", y = "ECDF", title = ttl) +
    theme_bw(base_size = BASE) +
    theme(
      aspect.ratio = 1, legend.position = c(0.02, 0.99), legend.justification = c(0, 1),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.key = element_rect(fill = "transparent", colour = NA)
    )
}
ecdf_B <- rbindlist(list(
  data.table(size_mb = sim_truth_sz, series = L_st),
  data.table(size_mb = sim_sz_sim, series = L_sim),
  data.table(size_mb = sim_sz_molb, series = L_molb)
))
ecdf_C <- rbindlist(list(
  data.table(size_mb = molb_truth_sz, series = L_mt),
  data.table(size_mb = skim_sz_sim, series = L_sim),
  data.table(size_mb = skim_sz_molb, series = L_molb)
))
p_B <- mk_ecdf(ecdf_B, "Simulation:\nnnil calls vs latent ancestry", c(L_st, L_sim, L_molb))
p_C <- mk_ecdf(ecdf_C, "Real skim:\nnnil calls vs MolBreeding truth", c(L_mt, L_molb, L_sim))

qs <- ppoints(200)
qq <- data.table(sim = quantile(skim_sz_sim, qs), molb = quantile(skim_sz_molb, qs))
p_D <- ggplot(qq, aes(molb, sim)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50") +
  geom_point(size = 1.3, colour = "grey20") +
  scale_x_log10(limits = xlim_mb, breaks = MB_BREAKS, oob = scales::oob_keep) +
  scale_y_log10(limits = xlim_mb, breaks = MB_BREAKS, oob = scales::oob_keep) +
  labs(
    x = "Introgression size(Mb)\nfrom molbio calibrated calls",
    y = "Introgression size(Mb)\nfrom sim calibrated calls",
    title = "QQ: sim- vs MolBreeding-\ncalibrated nnil (real skim)"
  ) +
  theme_bw(base_size = BASE) +
  theme(aspect.ratio = 1)

fig <- (p_A | p_B) / (p_C | p_D) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 25, face = "bold"),
    plot.tag.location = "plot", plot.tag.position = "topleft"
  )
ggsave(file.path(OUT, "fig_zeal_molbreeding_skim_calibration.png"), fig, width = 13, height = 13, dpi = 150)
# Figure caption (for the notebook): panel A open circles = each calibration's marker-mismatch
# optimum; dashed vertical = taxon nir (mean REF+0.5*HET over all taxa).
log_info("[fig] wrote png | caption note: A circle = mismatch-opt, dashed vertical = taxon nir (%.2f)", nir_taxon)
