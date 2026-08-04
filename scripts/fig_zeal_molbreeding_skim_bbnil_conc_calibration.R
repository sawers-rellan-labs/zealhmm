#!/usr/bin/env Rscript
# ZEAL molbreeding-skim bbnil CONC calibration figure -- the emission-axis analog of
# fig_zeal_molbreeding_skim_bbnil_calibration.R (which swept rrate). Motivation: rrate is
# the geometric self-transition, i.e. the fragment-LENGTH prior, so calibrating it against a
# fragment-size KS is near-circular, and it is a physical rate the genetic map already pins.
# The genuine skim-coverage mis-specification is in the EMISSION: how noisy, shallow, ref-biased
# read counts map onto REF/HET/ALT. For bbnil (count/BetaBinomial emission) that knob is the
# concentration `conc` (overdispersion) -- the count-emission analog of nnil's `nir`.
#
# So here: rrate is FIXED at the map-derived value (sim = SNP50K map_r; molb = R_MOLB), means are
# EM-FIT (fit_means = TRUE, the count-emission extension), and `conc` is the swept knob. With the
# fragment-length prior fixed, the fragment-size KS becomes a legitimate criterion (only the
# emission moves). fit_means = TRUE => each conc is a distinct EM fit (no redundancy), so this
# routes through call_ancestry per conc (caller_grid forbids fit_means).
#
# Truths are the SAME as the rrate figure:
#   SIM  = sim skim counts vs the simulated latent ancestry (BC2S3).
#   MolB = real skim counts vs nnil-on-molbreeding-hardcalls (the independent truth).
# make_cal_fig(crit): "mismatch", "ks" (introgression-size KS vs the genotype-derived ancestry).
#
# Compute CACHED; NNIL_ZEALBBC_RECOMPUTE=1 forces a rebuild.
#   Rscript scripts/fig_zeal_molbreeding_skim_bbnil_conc_calibration.R

suppressMessages({
  library(devtools)
  library(data.table)
  library(here)
  library(ggplot2)
  library(patchwork)
  library(ggtext)
})
devtools::load_all(path.expand("~/repos/nilhmm"), quiet = TRUE) # call_ancestry (feat branch)
source(here::here("R/metrics.R"))
source(here::here("scripts/logging.R"))

OUT <- here::here("results/sim/zeal_nil")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CONC_GRID <- c(2, 5, 10, 20, 50, 100, 200, 500) # BetaBinomial concentration (low = overdispersed)
ERR <- 0.01 # bbnil per-read error (fixed; initialises the EM means)
R_MOLB <- 1.67e-3 # MAP-derived rate: molb truth caller AND the molb skim test decode
TRUTH_CFG <- data.frame(nir = 0.7, germ = 1e-3, gert = 1e-4, p = 0.9, mr = 0)
KEEP <- c("name", "chr", "start_bp", "end_bp", "state")
BASE <- 20
ANN <- BASE * 0.8 / .pt
CACHE <- file.path(OUT, "fig_zeal_molbreeding_skim_bbnil_conc_calibration_cache.rds")

log_info(
  "[bbc] emission means: EM-fit (fit_means=TRUE); rrate FIXED at map value; sweeping conc %s",
  paste(CONC_GRID, collapse = ",")
)

if (file.exists(CACHE) && Sys.getenv("NNIL_ZEALBBC_RECOMPUTE") == "") {
  log_info("[bbc] reusing cached compute (%s)", basename(CACHE))
  D <- readRDS(CACHE)
} else {
  seg_at <- function(seg, v) seg[abs(conc - v) < 1e-9]
  keyc <- function(v) sprintf("%g", v)
  # one bbnil decode of the whole cohort at a fixed conc + fixed map rrate (fit_means = TRUE).
  # caller_grid forbids fit_means (emission is data-coupled), so route through call_ancestry.
  bb_conc <- function(dat, cc, rr, threads) {
    s <- as.data.table(as.data.frame(call_ancestry(as.data.frame(dat),
      caller = "bbnil", fit_means = TRUE, conc = cc, err = ERR,
      rrate = rr, design = "BC2S3", parallel = threads > 1L, threads = threads
    )))
    s$conc <- cc
    s
  }

  # ---- SIM leg: bbnil on sim skim COUNTS, sweep conc at map rrate; truth = latent ancestry ----
  sim <- readRDS(file.path(OUT, "zeal_nil_bc2s3_full.rds"))
  N_SIM <- as.integer(Sys.getenv("NNIL_ZEALBBC_NSIM", "500"))
  set.seed(1)
  sub <- sort(sample(seq_along(sim$names), min(N_SIM, length(sim$names))))
  sub_names <- sim$names[sub]
  M <- nrow(sim$grid)
  sim_long <- data.table(
    name = rep(sub_names, each = M),
    chr = rep(as.integer(sim$grid$chr), length(sub)),
    pos = rep(as.integer(sim$grid$pos), length(sub)),
    n_ref = as.integer(sim$n_ref[, sub]),
    n_alt = as.integer(sim$n_alt[, sub])
  )
  sim_truth <- as.data.table(sim$truth)[name %in% sub_names]
  sim_grid <- data.table(chr = as.integer(sim$grid$chr), pos = as.integer(sim$grid$pos))
  sim_traster <- rasterize_named(sim_truth, sim_grid)
  sim_tblocks <- .donor_blocks(sim_truth)
  R_SKIM <- sim$map_r # SNP50K map-derived rate: the FIXED sim rrate (test decode + gtrue ref)
  threads <- min(parallel::detectCores() - 2L, 8L)
  log_info(
    "[bbc] SIM: %d/%d lines x %d markers; fixed rrate=%.3e err=%g; conc grid %s",
    length(sub), length(sim$names), M, R_SKIM, ERR, paste(CONC_GRID, collapse = ",")
  )
  t0 <- Sys.time()
  sim_seg <- rbindlist(lapply(CONC_GRID, function(cc) {
    s <- bb_conc(sim_long, cc, R_SKIM, threads)
    log_info("[bbc] sim conc=%g elapsed %.0fs", cc, as.numeric(difftime(Sys.time(), t0, units = "secs")))
    s
  }))
  sim_mm <- sapply(CONC_GRID, function(v) 1 - marker_dsc(seg_at(sim_seg, v)[, ..KEEP], sim_truth, sim_grid, truth_raster = sim_traster)$accuracy)
  sim_dsc <- sapply(CONC_GRID, function(v) donor_fragment_dsc(seg_at(sim_seg, v)[, ..KEEP], sim_truth, truth_blocks = sim_tblocks)$dsc)

  # ---- MolB truth (nnil on molbreeding hard calls); skim test: bbnil on COUNTS at R_MOLB ----
  pair <- fread(here::here("data/zeal/correspondence/calibration_pairing.csv"))[in_calibration == TRUE]
  hc <- fread(here::here("data/zeal/molbreeding/molbreeding_hardcalls_wsfilt.tsv"))[name %in% pair$truth_sample]
  hc[is.na(g), g := 3L]
  molb_grid <- unique(hc[, .(chr, pos)])[order(chr, pos)]
  molb_truth <- as.data.table(caller_grid(hc[, .(name, chr, pos, g)],
    caller = "nnil", emission_grid = TRUTH_CFG, rrate = R_MOLB, design = "BC2S3", threads = 4L
  ))[, ..KEEP]
  molb_truth[, name := pair$test_sample[match(name, pair$truth_sample)]]
  molb_traster <- rasterize_named(molb_truth, molb_grid)
  molb_tblocks <- .donor_blocks(molb_truth)

  skim <- rbindlist(lapply(pair$test_sample, function(s) {
    cf <- fread(here::here("data/zeal/skim/counts", paste0(s, ".tsv")),
      header = FALSE,
      col.names = c("contig", "pos", "rb", "rc", "ab", "ac")
    )
    data.table(name = s, chr = as.integer(sub("chr", "", cf$contig)), pos = cf$pos, n_ref = cf$rc, n_alt = cf$ac)
  }))
  setorder(skim, name, chr, pos)
  skim_seg <- rbindlist(lapply(CONC_GRID, function(cc) bb_conc(skim[, .(name, chr, pos, n_ref, n_alt)], cc, R_MOLB, 4L)))
  skim_mm <- sapply(CONC_GRID, function(v) 1 - marker_dsc(seg_at(skim_seg, v)[, ..KEEP], molb_truth, molb_grid, truth_raster = molb_traster)$accuracy)
  skim_dsc <- sapply(CONC_GRID, function(v) donor_fragment_dsc(seg_at(skim_seg, v)[, ..KEEP], molb_truth, truth_blocks = molb_tblocks)$dsc)

  # KS reference: genotype-derived ancestry (nnil on g_true / hardcalls), never a bbnil run
  sim_gtrue_long <- data.table(
    name = rep(sub_names, each = M), chr = rep(as.integer(sim$grid$chr), length(sub)),
    pos = rep(as.integer(sim$grid$pos), length(sub)), g = as.integer(sim$g_true[, sub])
  )
  sim_gref_sz <- donor_block_sizes(as.data.table(caller_grid(sim_gtrue_long,
    caller = "nnil",
    emission_grid = TRUTH_CFG, rrate = R_SKIM, design = "BC2S3", threads = threads
  ))[, ..KEEP])
  molb_ref_sz <- donor_block_sizes(molb_truth)
  sim_ks <- sapply(CONC_GRID, function(v) fragment_size_ks(donor_block_sizes(seg_at(sim_seg, v)[, ..KEEP]), sim_gref_sz))
  skim_ks <- sapply(CONC_GRID, function(v) fragment_size_ks(donor_block_sizes(seg_at(skim_seg, v)[, ..KEEP]), molb_ref_sz))

  D <- list(
    sweep = rbind(
      data.table(truth = "sim", conc = CONC_GRID, mismatch = sim_mm, dsc = sim_dsc, ks = sim_ks),
      data.table(truth = "molb", conc = CONC_GRID, mismatch = skim_mm, dsc = skim_dsc, ks = skim_ks)
    ),
    sim_truth_sz = donor_block_sizes(sim_truth),
    sim_gref_sz = sim_gref_sz,
    molb_truth_sz = donor_block_sizes(molb_truth),
    sim_sz_by_conc = setNames(lapply(CONC_GRID, function(v) donor_block_sizes(seg_at(sim_seg, v)[, ..KEEP])), keyc(CONC_GRID)),
    skim_sz_by_conc = setNames(lapply(CONC_GRID, function(v) donor_block_sizes(seg_at(skim_seg, v)[, ..KEEP])), keyc(CONC_GRID))
  )
  saveRDS(D, CACHE)
  fwrite(D$sweep, file.path(OUT, "fig_zeal_molbreeding_skim_bbnil_conc_calibration_sweep.csv"))
}
list2env(D, environment())

# ============================== figure builder ==============================
col_sim <- "#0072B2"
col_molb <- "#D55E00"
MB_BREAKS <- c(0.1, 1, 10, 100)
keyc <- function(v) sprintf("%g", v)

CRIT <- list(
  mismatch = list(
    metric = "mismatch", opt = which.min, ylab = "marker mismatch rate",
    title = "bbnil conc calibration (fit means, r=map)", legA = c(0.5, 0.99), legAj = c(0.5, 1),
    sim_ref = "sim_truth_sz", st_lab = "sim latent ancestry\nBC2S3",
    out = "fig_zeal_molbreeding_skim_bbnil_conc_calibration.png"
  ),
  ks = list(
    metric = "ks", opt = which.min, ylab = "introgression-size KS distance (D)",
    title = "bbnil conc KS calibration (fit means, r=map)", legA = c(0.5, 0.99), legAj = c(0.5, 1),
    sim_ref = "sim_gref_sz", st_lab = "sim nnil-on-g_true\nancestry",
    out = "fig_zeal_molbreeding_skim_bbnil_conc_ks_calibration.png"
  )
)

make_cal_fig <- function(crit) {
  cf <- CRIT[[crit]]
  sw <- copy(sweep)[, val := get(cf$metric)]
  cc_sim <- sw[truth == "sim"]$conc[cf$opt(sw[truth == "sim"]$val)]
  cc_molb <- sw[truth == "molb"]$conc[cf$opt(sw[truth == "molb"]$val)]
  log_info("[bbc:%s] operating points: sim conc=%g | molb conc=%g", crit, cc_sim, cc_molb)

  sim_sz_sim <- sim_sz_by_conc[[keyc(cc_sim)]]
  sim_sz_molb <- sim_sz_by_conc[[keyc(cc_molb)]]
  skim_sz_sim <- skim_sz_by_conc[[keyc(cc_sim)]]
  skim_sz_molb <- skim_sz_by_conc[[keyc(cc_molb)]]

  L_sim <- sprintf("sim calibrated \nbbnil (conc=%g)", cc_sim)
  L_molb <- sprintf("molb calibrated \nbbnil (conc=%g)", cc_molb)
  st_sz <- get(cf$sim_ref)
  L_st <- cf$st_lab
  L_mt <- "molb calls"
  pal <- c(setNames(c(col_sim, col_molb), c(L_sim, L_molb)), setNames(c("black", "black"), c(L_st, L_mt)))
  lty <- c(setNames(c("solid", "solid"), c(L_sim, L_molb)), setNames(c("dotted", "dotted"), c(L_st, L_mt)))
  allsz <- c(st_sz, sim_sz_sim, sim_sz_molb, molb_truth_sz, skim_sz_sim, skim_sz_molb)
  xlim_mb <- c(0.1, max(allsz[is.finite(allsz) & allsz > 0]))

  optrow <- sw[, .SD[cf$opt(val)], by = truth]
  p_A <- ggplot(sw, aes(conc, val, colour = truth)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.4) +
    geom_point(data = optrow, size = 3.6, shape = 21, fill = "white", stroke = 1.1) +
    scale_x_log10() +
    expand_limits(y = 0) +
    scale_colour_manual(
      values = c(sim = col_sim, molb = col_molb),
      labels = c(sim = "sim skim vs sim ancestry", molb = "real skim vs molb ancestry"), name = "calibration"
    ) +
    labs(x = "bbnil BetaBinomial conc (overdispersion)", y = cf$ylab, title = cf$title) +
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
    data.table(size_mb = st_sz, series = L_st),
    data.table(size_mb = sim_sz_sim, series = L_sim),
    data.table(size_mb = sim_sz_molb, series = L_molb)
  ))
  ecdf_C <- rbindlist(list(
    data.table(size_mb = molb_truth_sz, series = L_mt),
    data.table(size_mb = skim_sz_sim, series = L_sim),
    data.table(size_mb = skim_sz_molb, series = L_molb)
  ))
  p_B <- mk_ecdf(ecdf_B, "Simulation:\nbbnil calls vs latent ancestry", c(L_st, L_sim, L_molb))
  p_C <- mk_ecdf(ecdf_C, "Real skim:\nbbnil calls vs MolBreeding truth", c(L_mt, L_molb, L_sim))

  ks_C <- suppressWarnings(ks.test(skim_sz_sim, skim_sz_molb))
  ks_txt <- if (ks_C$p.value >= 0.01) sprintf("%.2f", ks_C$p.value) else sprintf("%.0e", ks_C$p.value)
  log_info("[bbc:%s] panel C KS (sim-cal vs molb-cal): D=%.3f p=%.3g", crit, ks_C$statistic, ks_C$p.value)
  ks_rich <- sprintf(
    "calibration<br><span style='color:%s'>molb</span> vs <span style='color:%s'>sim</span><br>*p* = %s",
    col_molb, col_sim, ks_txt
  )
  p_C <- p_C + geom_richtext(
    data = data.frame(x = 20, y = 0.25, label = ks_rich),
    aes(x, y, label = label), inherit.aes = FALSE,
    hjust = 0, vjust = 0.5, size = ANN, lineheight = 0.9, colour = "grey20",
    fill = "transparent", label.color = NA, label.r = grid::unit(0, "pt"), label.padding = grid::unit(2, "pt")
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
      title = "QQ: sim- vs MolBreeding-\ncalibrated bbnil (real skim)"
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
  log_info("[bbc:%s] wrote %s", crit, cf$out)
}

make_cal_fig("mismatch")
make_cal_fig("ks")
