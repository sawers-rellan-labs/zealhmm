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
# ZEAL has TWO calibration criteria; make_cal_fig() emits ONE figure per criterion:
#   mismatch  nir* = argmin marker mismatch (Holland's criterion) -> ..._calibration.png
#   fragment  nir* = argmax donor-fragment DSC                    -> ..._fragment_calibration.png
#
# Panels (2x2), per figure -- the SIM operating point and the MolB operating point of the
# chosen criterion are carried through B/C/D:
#   A  criterion vs nir, nnil scored against the SIM truth (blue) and the MolB truth (orange);
#      open circles = each truth's operating point;
#      dashed vertical = taxon nir ~0.69 (mean REF+0.5*HET over all taxa).
#   B  SIM: introgression-size ECDF, nnil at the two operating points vs the sim latent truth.
#   C  Real skim: introgression-size ECDF, nnil at the two operating points vs the MolB truth;
#      ggtext KS p-value (sim- vs molb-calibrated skim sizes).
#   D  QQ of real-skim introgression size, sim- (y) vs molb-calibrated (x).
#
# Compute is CACHED and criterion-agnostic; NNIL_ZEALCAL_RECOMPUTE=1 forces a rebuild.
#   Rscript scripts/fig_zeal_molbreeding_skim_calibration.R

suppressMessages({
  library(devtools)
  library(data.table)
  library(here)
  library(ggplot2)
  library(patchwork)
  library(ggtext)
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
ANN <- BASE * 0.8 / .pt # match the legend category-label size (theme legend.text = rel(0.8) of BASE)
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
  sim_tblocks <- .donor_blocks(sim_truth) # merge sim truth blocks ONCE (for the DSC criterion)
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
  sim_dsc <- sapply(NIR_GRID, function(v) donor_fragment_dice(seg_at(sim_seg, v)[, ..KEEP], sim_truth, truth_blocks = sim_tblocks)$dice)

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

  # Criterion-agnostic cache: store BOTH calibration curves (marker mismatch AND donor-
  # fragment DSC) for each truth, and the introgression-size distribution at EVERY nir grid
  # point (sim + real skim). make_cal_fig() then picks the operating point for either
  # criterion and looks up the matching size distributions -- no re-decoding per figure.
  keyv <- function(v) sprintf("%.2f", v)
  D <- list(
    sweep = rbind(
      data.table(truth = "sim", nir = NIR_GRID, mismatch = sim_mm, dsc = sim_dsc),
      data.table(truth = "molb", nir = NIR_GRID, mismatch = skim_mm, dsc = skim_dsc)
    ),
    sim_truth_sz = donor_block_sizes(sim_truth),
    molb_truth_sz = donor_block_sizes(molb_truth),
    sim_sz_by_nir = setNames(lapply(NIR_GRID, function(v) donor_block_sizes(seg_at(sim_seg, v)[, ..KEEP])), keyv(NIR_GRID)),
    skim_sz_by_nir = setNames(lapply(NIR_GRID, function(v) donor_block_sizes(seg_at(skim_seg, v)[, ..KEEP])), keyv(NIR_GRID))
  )
  saveRDS(D, CACHE)
  fwrite(D$sweep, file.path(OUT, "fig_zeal_molbreeding_skim_calibration_sweep.csv"))
}
list2env(D, environment())

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

# ============================== figure builder ==============================
# colour is CONSISTENT across all panels: blue = sim (calibration source / sim-calibrated
# caller), orange = MolBreeding (source / MolBreeding-calibrated caller), black dotted = truth.
# make_cal_fig(crit) builds the same 4-panel figure for EITHER calibration criterion:
#   "mismatch" = argmin marker mismatch (Holland's criterion)
#   "fragment" = argmax donor-fragment DSC
# Interior legend / KS-annotation placement is set per criterion (the mismatch curve is a
# valley, the DSC curve rises to the right, so the empty corner differs).
col_sim <- "#0072B2"
col_molb <- "#D55E00"
MB_BREAKS <- c(0.1, 1, 10, 100)
keyv <- function(v) sprintf("%.2f", v)

CRIT <- list(
  mismatch = list(
    metric = "mismatch", opt = which.min, ylab = "marker mismatch rate",
    title = "nnil nir mismatch calibration", legA = c(0.02, 0.99), legAj = c(0, 1),
    taxon_y = 0.30, ks_x = 20, ks_y = 0.25, ks_hjust = 0,
    out = "fig_zeal_molbreeding_skim_calibration.png"
  ),
  fragment = list(
    metric = "dsc", opt = which.max, ylab = "donor fragment DSC",
    title = "nnil nir fragment calibration", legA = c(0.02, 0.99), legAj = c(0, 1),
    taxon_y = 0.30, ks_x = 20, ks_y = 0.25, ks_hjust = 0,
    out = "fig_zeal_molbreeding_skim_fragment_calibration.png"
  )
)

make_cal_fig <- function(crit) {
  cf <- CRIT[[crit]]
  sw <- copy(sweep)[, val := get(cf$metric)]
  nir_sim <- sw[truth == "sim"]$nir[cf$opt(sw[truth == "sim"]$val)]
  nir_molb <- sw[truth == "molb"]$nir[cf$opt(sw[truth == "molb"]$val)]
  log_info("[fig:%s] operating points: sim nir=%.2f | molb nir=%.2f", crit, nir_sim, nir_molb)

  sim_sz_sim <- sim_sz_by_nir[[keyv(nir_sim)]]
  sim_sz_molb <- sim_sz_by_nir[[keyv(nir_molb)]]
  skim_sz_sim <- skim_sz_by_nir[[keyv(nir_sim)]]
  skim_sz_molb <- skim_sz_by_nir[[keyv(nir_molb)]]

  L_sim <- sprintf("sim calibrated \nnnil (nir=%.2f)", nir_sim)
  L_molb <- sprintf("molb calibrated \nnnil (nir=%.2f)", nir_molb)
  L_st <- "sim latent ancestry\nBC2S3"
  L_mt <- "molb calls"
  pal <- c(setNames(c(col_sim, col_molb), c(L_sim, L_molb)), setNames(c("black", "black"), c(L_st, L_mt)))
  lty <- c(setNames(c("solid", "solid"), c(L_sim, L_molb)), setNames(c("dotted", "dotted"), c(L_st, L_mt)))
  allsz <- c(sim_truth_sz, sim_sz_sim, sim_sz_molb, molb_truth_sz, skim_sz_sim, skim_sz_molb)
  xlim_mb <- c(0.1, max(allsz[is.finite(allsz) & allsz > 0]))

  # A: sim (blue) + MolBreeding (orange) calibration curves; open circle = operating point
  optrow <- sw[, .SD[cf$opt(val)], by = truth]
  p_A <- ggplot(sw, aes(nir, val, colour = truth)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.4) +
    geom_point(data = optrow, size = 3.6, shape = 21, fill = "white", stroke = 1.1) +
    geom_vline(xintercept = nir_taxon, linetype = "dashed", colour = "black", linewidth = 0.5, alpha = 0.45) +
    annotate("text",
      x = nir_taxon, y = cf$taxon_y, label = "taxon",
      angle = 90, hjust = 0.5, vjust = 1.4, size = BASE * 0.24, colour = "black", alpha = 0.45
    ) +
    scale_colour_manual(
      values = c(sim = col_sim, molb = col_molb),
      labels = c(sim = "sim skim vs sim ancestry", molb = "real skim vs molb ancestry"), name = "calibration"
    ) +
    labs(x = expression("non-informative rate " * italic(nir)), y = cf$ylab, title = cf$title) +
    theme_bw(base_size = BASE) +
    theme(
      aspect.ratio = 1, legend.position = cf$legA, legend.justification = cf$legAj,
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.key = element_rect(fill = "transparent", colour = NA)
    )

  mk_ecdf <- function(d, ttl, lvls) {
    d <- copy(d)[, series := factor(series, levels = lvls)]
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

  # KS on panel C: sim-calibrated vs molb-calibrated skim introgression sizes (same 14 NILs)
  ks_C <- suppressWarnings(ks.test(skim_sz_sim, skim_sz_molb))
  ks_txt <- if (ks_C$p.value >= 0.01) sprintf("%.2f", ks_C$p.value) else sprintf("%.0e", ks_C$p.value)
  log_info("[fig:%s] panel C KS (sim-cal vs molb-cal nnil): D=%.3f p=%.3g", crit, ks_C$statistic, ks_C$p.value)
  ks_rich <- sprintf(
    "calibration<br><span style='color:%s'>molb</span> vs <span style='color:%s'>sim</span><br>*p* = %s",
    col_molb, col_sim, ks_txt
  )
  p_C <- p_C + geom_richtext(
    data = data.frame(x = cf$ks_x, y = cf$ks_y, label = ks_rich),
    aes(x, y, label = label), inherit.aes = FALSE,
    hjust = cf$ks_hjust, vjust = 0.5, size = ANN, lineheight = 0.9, colour = "grey20",
    fill = "transparent", label.color = NA, label.r = grid::unit(0, "pt"),
    label.padding = grid::unit(2, "pt")
  )

  qs <- ppoints(200)
  qq <- data.table(sim = quantile(skim_sz_sim, qs), molb = quantile(skim_sz_molb, qs))
  p_D <- ggplot(qq, aes(molb, sim)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey50") +
    geom_point(size = 1.3, colour = col_sim) +
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
  ggsave(file.path(OUT, cf$out), fig, width = 13, height = 13, dpi = 150)
  log_info("[fig:%s] wrote %s (taxon nir=%.2f)", crit, cf$out, nir_taxon)
}

make_cal_fig("mismatch")
make_cal_fig("fragment")
