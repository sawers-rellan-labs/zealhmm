#!/usr/bin/env Rscript
# ZEAL molbreeding-skim rtiger RIGIDITY calibration figure -- the rtiger analog of
# fig_zeal_molbreeding_skim_calibration.R (nnil nir). Baseline = the mismatch figure.
#
# rtiger consumes read counts (n_ref/n_alt) directly (NO hard-calling) and its only
# calibrated knob here is the RIGIDITY (minimum run length); caller_sweep(refit="none")
# fits the count/BetaBinomial emission ONCE and decodes per rigidity. Truths are the
# SAME as the nnil figure:
#   SIM  = sim skim counts vs the simulated latent ancestry (BC2S3).
#   MolB = real skim counts vs nnil-on-molbreeding-hardcalls (the independent truth).
# Criterion = marker mismatch (argmin), scored identically to the nnil figure.
#
# Panels (2x2): A mismatch vs rigidity (log x), sim (blue) + MolB (orange), open circle
# = argmin; B sim introgression-size ECDF at the two operating points vs latent ancestry;
# C real-skim ECDF vs MolB truth + ggtext KS (sim- vs molb-calibrated); D QQ (sim y / molb x).
# There is NO taxon-nir reference line: rigidity is a segmentation length, not nir.
#
# Compute CACHED; NNIL_ZEALRTIG_RECOMPUTE=1 forces a rebuild.
#   Rscript scripts/fig_zeal_molbreeding_skim_rtiger_calibration.R

suppressMessages({
  library(devtools)
  library(data.table)
  library(here)
  library(ggplot2)
  library(patchwork)
  library(ggtext)
})
devtools::load_all(path.expand("~/repos/nilhmm"), quiet = TRUE) # caller_sweep, caller_grid (feat branch)
source(here::here("R/metrics.R"))
source(here::here("scripts/logging.R"))

OUT <- here::here("results/sim/zeal_nil")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
RIG_GRID <- c(2L, 5L, 10L, 20L, 40L, 60L, 80L) # min run length; max 80 < 2*min-covered floor (sim chr9 = 193)
R_MOLB <- 1.67e-3
TRUTH_CFG <- data.frame(nir = 0.7, germ = 1e-3, gert = 1e-4, p = 0.9, mr = 0) # MolB truth caller (nnil on hardcalls)
KEEP <- c("name", "chr", "start_bp", "end_bp", "state")
BASE <- 20
ANN <- BASE * 0.8 / .pt
CACHE <- file.path(OUT, "fig_zeal_molbreeding_skim_rtiger_calibration_cache.rds")

if (file.exists(CACHE) && Sys.getenv("NNIL_ZEALRTIG_RECOMPUTE") == "") {
  log_info("[rtig] reusing cached compute (%s)", basename(CACHE))
  D <- readRDS(CACHE)
} else {
  ml_g <- function(nr, na) {
    g <- nilHMM::call_gt(as.integer(nr), as.integer(na), prior = "flat", error = 0.01, return = "call")
    g[is.na(g)] <- 3L
    as.integer(g)
  }
  seg_at <- function(seg, v) seg[rigidity == v]
  keyv <- function(v) sprintf("%d", as.integer(v))

  # ---- SIM leg: rtiger on sim skim COUNTS; truth = latent ancestry ----
  sim <- readRDS(file.path(OUT, "zeal_nil_bc2s3_full.rds"))
  N_SIM <- as.integer(Sys.getenv("NNIL_ZEALRTIG_NSIM", "500"))
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
  threads <- min(parallel::detectCores() - 2L, 8L)
  log_info("[rtig] SIM: %d/%d lines x %d markers; rigidity grid %s", length(sub), length(sim$names), M, paste(RIG_GRID, collapse = ","))
  t0 <- Sys.time()
  sim_seg <- as.data.table(caller_sweep(sim_long, caller = "rtiger", values = RIG_GRID, refit = "none", design = "BC2S3", threads = threads))
  log_info("[rtig] SIM rtiger sweep done in %.0fs", as.numeric(difftime(Sys.time(), t0, units = "secs")))
  # full metric suite per rigidity (all kept in the sweep so any criterion is calibratable)
  metrics_leg <- function(seg, truth, grid, tr, tb) {
    rbindlist(lapply(RIG_GRID, function(v) {
      s <- seg_at(seg, v)[, ..KEEP]
      md <- marker_dice(s, truth, grid, truth_raster = tr)
      dr <- md$per_class[class == "donor(>0)"]
      ff <- donor_fragment_dice(s, truth, truth_blocks = tb)
      data.table(
        rigidity = v, mismatch = 1 - md$accuracy, mk_recall = dr$recall[1],
        dsc = ff$dice, fdr = ff$fdr, frag_recall = ff$recall,
        n_called = ff$n_called, n_truth = ff$n_truth, nbreak = breakpoint_count(s)
      )
    }))
  }
  sim_m <- metrics_leg(sim_seg, sim_truth, sim_grid, sim_traster, sim_tblocks)

  # ---- MolB truth (nnil on molbreeding hard calls); skim test: rtiger on COUNTS ----
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
  skim_seg <- as.data.table(caller_sweep(skim[, .(name, chr, pos, n_ref, n_alt)],
    caller = "rtiger", values = RIG_GRID, refit = "none", design = "BC2S3", threads = 4L
  ))
  skim_m <- metrics_leg(skim_seg, molb_truth, molb_grid, molb_traster, molb_tblocks)

  D <- list(
    sweep = rbind(cbind(truth = "sim", sim_m), cbind(truth = "molb", skim_m)),
    sim_truth_sz = donor_block_sizes(sim_truth),
    molb_truth_sz = donor_block_sizes(molb_truth),
    sim_sz_by_rig = setNames(lapply(RIG_GRID, function(v) donor_block_sizes(seg_at(sim_seg, v)[, ..KEEP])), keyv(RIG_GRID)),
    skim_sz_by_rig = setNames(lapply(RIG_GRID, function(v) donor_block_sizes(seg_at(skim_seg, v)[, ..KEEP])), keyv(RIG_GRID))
  )
  saveRDS(D, CACHE)
  fwrite(D$sweep, file.path(OUT, "fig_zeal_molbreeding_skim_rtiger_calibration_sweep.csv"))
}
list2env(D, environment())

# ============================== figure builder ==============================
# make_cal_fig(crit) builds the 4-panel figure calibrating rtiger rigidity on `crit`
# (all metrics live in the sweep). "mismatch" = argmin marker mismatch (baseline);
# "fdr" = argmin donor-fragment FDR. Panels B/C/D use each leg's chosen rigidity.
col_sim <- "#0072B2"
col_molb <- "#D55E00"
MB_BREAKS <- c(0.1, 1, 10, 100)
keyv <- function(v) sprintf("%d", as.integer(v))

CRIT <- list(
  mismatch = list(
    metric = "mismatch", opt = which.min, ylab = "marker mismatch rate",
    title = "rtiger rigidity calibration", legA = c(0.02, 0.02), legAj = c(0, 0),
    out = "fig_zeal_molbreeding_skim_rtiger_calibration.png"
  ),
  fdr = list(
    metric = "fdr", opt = which.min, ylab = "donor fragment FDR",
    title = "rtiger rigidity FDR calibration", legA = c(0.98, 0.98), legAj = c(1, 1),
    out = "fig_zeal_molbreeding_skim_rtiger_fdr_calibration.png"
  )
)

make_cal_fig <- function(crit) {
  cf <- CRIT[[crit]]
  sw <- copy(sweep)[, val := get(cf$metric)]
  rig_sim <- sw[truth == "sim"]$rigidity[cf$opt(sw[truth == "sim"]$val)]
  rig_molb <- sw[truth == "molb"]$rigidity[cf$opt(sw[truth == "molb"]$val)]
  log_info("[rtig:%s] operating points: sim R=%d | molb R=%d", crit, rig_sim, rig_molb)

  sim_sz_sim <- sim_sz_by_rig[[keyv(rig_sim)]]
  sim_sz_molb <- sim_sz_by_rig[[keyv(rig_molb)]]
  skim_sz_sim <- skim_sz_by_rig[[keyv(rig_sim)]]
  skim_sz_molb <- skim_sz_by_rig[[keyv(rig_molb)]]

  L_sim <- sprintf("sim calibrated \nrtiger (R=%d)", rig_sim)
  L_molb <- sprintf("molb calibrated \nrtiger (R=%d)", rig_molb)
  L_st <- "sim latent ancestry\nBC2S3"
  L_mt <- "molb calls"
  pal <- c(setNames(c(col_sim, col_molb), c(L_sim, L_molb)), setNames(c("black", "black"), c(L_st, L_mt)))
  lty <- c(setNames(c("solid", "solid"), c(L_sim, L_molb)), setNames(c("dotted", "dotted"), c(L_st, L_mt)))
  allsz <- c(sim_truth_sz, sim_sz_sim, sim_sz_molb, molb_truth_sz, skim_sz_sim, skim_sz_molb)
  xlim_mb <- c(0.1, max(allsz[is.finite(allsz) & allsz > 0]))

  # A: criterion vs rigidity (log x), y from 0; open circle = operating point
  optrow <- sw[, .SD[cf$opt(val)], by = truth]
  p_A <- ggplot(sw, aes(rigidity, val, colour = truth)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.4) +
    geom_point(data = optrow, size = 3.6, shape = 21, fill = "white", stroke = 1.1) +
    scale_x_log10() +
    expand_limits(y = 0) +
    scale_colour_manual(
      values = c(sim = col_sim, molb = col_molb),
      labels = c(sim = "sim skim vs sim ancestry", molb = "real skim vs molb ancestry"), name = "calibration"
    ) +
    labs(x = "rtiger rigidity (min run length)", y = cf$ylab, title = cf$title) +
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
  p_B <- mk_ecdf(ecdf_B, "Simulation:\nrtiger calls vs latent ancestry", c(L_st, L_sim, L_molb))
  p_C <- mk_ecdf(ecdf_C, "Real skim:\nrtiger calls vs MolBreeding truth", c(L_mt, L_molb, L_sim))

  ks_C <- suppressWarnings(ks.test(skim_sz_sim, skim_sz_molb))
  ks_txt <- if (ks_C$p.value >= 0.01) sprintf("%.2f", ks_C$p.value) else sprintf("%.0e", ks_C$p.value)
  log_info("[rtig:%s] panel C KS (sim-cal vs molb-cal): D=%.3f p=%.3g", crit, ks_C$statistic, ks_C$p.value)
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
      title = "QQ: sim- vs MolBreeding-\ncalibrated rtiger (real skim)"
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
  log_info("[rtig:%s] wrote %s", crit, cf$out)
}

make_cal_fig("mismatch")
make_cal_fig("fdr")
