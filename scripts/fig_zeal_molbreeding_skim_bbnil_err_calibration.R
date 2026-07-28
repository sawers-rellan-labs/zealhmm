#!/usr/bin/env Rscript
# ZEAL molbreeding-skim bbnil ERR calibration figure -- the a-priori (Jim-faithful) emission
# calibration, the counterpart to the conc figure. With fit_means = FALSE, the count-emission
# state means are FIXED at theta = c(err, 0.5, 1 - err), so `err` IS the emission: it places the
# REF/ALT means by hand rather than fitting them. This is the "Holland does not fit means, he adds
# a priori parameters" path. rrate is FIXED at the map value (sim = map_r, molb = R_MOLB), so the
# fragment-length prior is pinned and only the emission moves; conc held at its default 20.
#
# Caveat: err is SYMMETRIC (shrinks REF and ALT toward 0.5 together), so it cannot represent
# reference bias (asymmetric REF/ALT means) -- the failure mode fit_means addresses. This figure
# tests whether a tuned a-priori err recovers the fit_means gain. fit_means = FALSE => the fast
# batched caller_grid C++ path (no EM), so this is cheap.
#
# Truths are the SAME as the conc/rrate figures:
#   SIM  = sim skim counts vs the simulated latent ancestry (BC2S3).
#   MolB = real skim counts vs nnil-on-molbreeding-hardcalls (the independent truth).
# make_cal_fig(crit): "mismatch", "ks" (introgression-size KS vs the genotype-derived ancestry).
#
# Compute CACHED; NNIL_ZEALBBE_RECOMPUTE=1 forces a rebuild.
#   Rscript scripts/fig_zeal_molbreeding_skim_bbnil_err_calibration.R

suppressMessages({
  library(devtools)
  library(data.table)
  library(here)
  library(ggplot2)
  library(patchwork)
  library(ggtext)
})
devtools::load_all(path.expand("~/repos/nilhmm"), quiet = TRUE) # caller_grid (feat branch)
source(here::here("R/metrics.R"))
source(here::here("scripts/logging.R"))

OUT <- here::here("results/sim/zeal_nil")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ERR_GRID <- c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2) # a-priori per-read error -> theta = c(err,0.5,1-err)
CONC_FIX <- 20 # BetaBinomial concentration held at its default
R_MOLB <- 1.67e-3 # MAP-derived rate: molb truth caller AND the molb skim test decode
TRUTH_CFG <- data.frame(nir = 0.7, germ = 1e-3, gert = 1e-4, p = 0.9, mr = 0)
KEEP <- c("name", "chr", "start_bp", "end_bp", "state")
BASE <- 20
ANN <- BASE * 0.8 / .pt
CACHE <- file.path(OUT, "fig_zeal_molbreeding_skim_bbnil_err_calibration_cache.rds")

log_info(
  "[bbe] emission means: FIXED c(err,0.5,1-err) (fit_means=FALSE, a-priori); rrate FIXED at map; sweeping err %s",
  paste(ERR_GRID, collapse = ",")
)

if (file.exists(CACHE) && Sys.getenv("NNIL_ZEALBBE_RECOMPUTE") == "") {
  log_info("[bbe] reusing cached compute (%s)", basename(CACHE))
  D <- readRDS(CACHE)
} else {
  seg_at <- function(seg, v) seg[abs(err - v) < v * 1e-6]
  keye <- function(v) sprintf("%g", v)
  # one bbnil decode of the whole cohort at a fixed err + fixed map rrate, means FIXED
  # (fit_means = FALSE) -> the batched caller_grid C++ path (no EM).
  bb_err <- function(dat, ee, rr, threads) {
    s <- as.data.table(caller_grid(dat,
      caller = "bbnil", emission_grid = data.frame(fit_means = FALSE),
      conc = CONC_FIX, err = ee, rrate = rr, design = "BC2S3", threads = threads
    ))
    s$err <- ee
    s
  }

  # ---- SIM leg: bbnil on sim skim COUNTS, sweep err at map rrate; truth = latent ancestry ----
  sim <- readRDS(file.path(OUT, "zeal_nil_bc2s3_full.rds"))
  N_SIM <- as.integer(Sys.getenv("NNIL_ZEALBBE_NSIM", "500"))
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
    "[bbe] SIM: %d/%d lines x %d markers; fixed rrate=%.3e conc=%g; err grid %s",
    length(sub), length(sim$names), M, R_SKIM, CONC_FIX, paste(ERR_GRID, collapse = ",")
  )
  t0 <- Sys.time()
  sim_seg <- rbindlist(lapply(ERR_GRID, function(ee) {
    s <- bb_err(sim_long, ee, R_SKIM, threads)
    log_info("[bbe] sim err=%g elapsed %.0fs", ee, as.numeric(difftime(Sys.time(), t0, units = "secs")))
    s
  }))
  sim_mm <- sapply(ERR_GRID, function(v) 1 - marker_dice(seg_at(sim_seg, v)[, ..KEEP], sim_truth, sim_grid, truth_raster = sim_traster)$accuracy)
  sim_dsc <- sapply(ERR_GRID, function(v) donor_fragment_dice(seg_at(sim_seg, v)[, ..KEEP], sim_truth, truth_blocks = sim_tblocks)$dice)

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
  skim_seg <- rbindlist(lapply(ERR_GRID, function(ee) bb_err(skim[, .(name, chr, pos, n_ref, n_alt)], ee, R_MOLB, 4L)))
  skim_mm <- sapply(ERR_GRID, function(v) 1 - marker_dice(seg_at(skim_seg, v)[, ..KEEP], molb_truth, molb_grid, truth_raster = molb_traster)$accuracy)
  skim_dsc <- sapply(ERR_GRID, function(v) donor_fragment_dice(seg_at(skim_seg, v)[, ..KEEP], molb_truth, truth_blocks = molb_tblocks)$dice)

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
  sim_ks <- sapply(ERR_GRID, function(v) fragment_size_ks(donor_block_sizes(seg_at(sim_seg, v)[, ..KEEP]), sim_gref_sz))
  skim_ks <- sapply(ERR_GRID, function(v) fragment_size_ks(donor_block_sizes(seg_at(skim_seg, v)[, ..KEEP]), molb_ref_sz))

  D <- list(
    sweep = rbind(
      data.table(truth = "sim", err = ERR_GRID, mismatch = sim_mm, dsc = sim_dsc, ks = sim_ks),
      data.table(truth = "molb", err = ERR_GRID, mismatch = skim_mm, dsc = skim_dsc, ks = skim_ks)
    ),
    sim_truth_sz = donor_block_sizes(sim_truth),
    sim_gref_sz = sim_gref_sz,
    molb_truth_sz = donor_block_sizes(molb_truth),
    sim_sz_by_err = setNames(lapply(ERR_GRID, function(v) donor_block_sizes(seg_at(sim_seg, v)[, ..KEEP])), keye(ERR_GRID)),
    skim_sz_by_err = setNames(lapply(ERR_GRID, function(v) donor_block_sizes(seg_at(skim_seg, v)[, ..KEEP])), keye(ERR_GRID))
  )
  saveRDS(D, CACHE)
  fwrite(D$sweep, file.path(OUT, "fig_zeal_molbreeding_skim_bbnil_err_calibration_sweep.csv"))
}
list2env(D, environment())

# ============================== figure builder ==============================
col_sim <- "#0072B2"
col_molb <- "#D55E00"
MB_BREAKS <- c(0.1, 1, 10, 100)
keye <- function(v) sprintf("%g", v)

CRIT <- list(
  mismatch = list(
    metric = "mismatch", opt = which.min, ylab = "marker mismatch rate",
    title = "bbnil err calibration (a-priori means, r=map)", legA = c(0.5, 0.99), legAj = c(0.5, 1),
    sim_ref = "sim_truth_sz", st_lab = "sim latent ancestry\nBC2S3",
    out = "fig_zeal_molbreeding_skim_bbnil_err_calibration.png"
  ),
  ks = list(
    metric = "ks", opt = which.min, ylab = "introgression-size KS distance (D)",
    title = "bbnil err KS calibration (a-priori means, r=map)", legA = c(0.5, 0.99), legAj = c(0.5, 1),
    sim_ref = "sim_gref_sz", st_lab = "sim nnil-on-g_true\nancestry",
    out = "fig_zeal_molbreeding_skim_bbnil_err_ks_calibration.png"
  )
)

make_cal_fig <- function(crit) {
  cf <- CRIT[[crit]]
  sw <- copy(sweep)[, val := get(cf$metric)]
  ee_sim <- sw[truth == "sim"]$err[cf$opt(sw[truth == "sim"]$val)]
  ee_molb <- sw[truth == "molb"]$err[cf$opt(sw[truth == "molb"]$val)]
  log_info("[bbe:%s] operating points: sim err=%g | molb err=%g", crit, ee_sim, ee_molb)

  sim_sz_sim <- sim_sz_by_err[[keye(ee_sim)]]
  sim_sz_molb <- sim_sz_by_err[[keye(ee_molb)]]
  skim_sz_sim <- skim_sz_by_err[[keye(ee_sim)]]
  skim_sz_molb <- skim_sz_by_err[[keye(ee_molb)]]

  L_sim <- sprintf("sim calibrated \nbbnil (err=%g)", ee_sim)
  L_molb <- sprintf("molb calibrated \nbbnil (err=%g)", ee_molb)
  st_sz <- get(cf$sim_ref)
  L_st <- cf$st_lab
  L_mt <- "molb calls"
  pal <- c(setNames(c(col_sim, col_molb), c(L_sim, L_molb)), setNames(c("black", "black"), c(L_st, L_mt)))
  lty <- c(setNames(c("solid", "solid"), c(L_sim, L_molb)), setNames(c("dotted", "dotted"), c(L_st, L_mt)))
  allsz <- c(st_sz, sim_sz_sim, sim_sz_molb, molb_truth_sz, skim_sz_sim, skim_sz_molb)
  xlim_mb <- c(0.1, max(allsz[is.finite(allsz) & allsz > 0]))

  optrow <- sw[, .SD[cf$opt(val)], by = truth]
  p_A <- ggplot(sw, aes(err, val, colour = truth)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.4) +
    geom_point(data = optrow, size = 3.6, shape = 21, fill = "white", stroke = 1.1) +
    scale_x_log10() +
    expand_limits(y = 0) +
    scale_colour_manual(
      values = c(sim = col_sim, molb = col_molb),
      labels = c(sim = "sim skim vs sim ancestry", molb = "real skim vs molb ancestry"), name = "calibration"
    ) +
    labs(x = "bbnil a-priori err (sets REF/ALT means)", y = cf$ylab, title = cf$title) +
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
  log_info("[bbe:%s] panel C KS (sim-cal vs molb-cal): D=%.3f p=%.3g", crit, ks_C$statistic, ks_C$p.value)
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
  log_info("[bbe:%s] wrote %s", crit, cf$out)
}

make_cal_fig("mismatch")
make_cal_fig("ks")
