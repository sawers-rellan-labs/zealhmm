#!/usr/bin/env Rscript
# ZEAL molbreeding-skim LB-Impute (lbimpute) RECOMBDIST calibration figure -- the
# lbimpute analog of fig_zeal_molbreeding_skim_calibration.R (nnil nir). Baseline = the
# mismatch figure.
#
# lbimpute is genetic-distance native: its transition decays over cM. We run it in cM
# (unit = "cm") on the SAME native v5 cM map the simulation was generated on
# (data/zeal/markers_snp50k_cm.tsv -- the map behind map_r), and calibrate its single
# knob, recombdist (the cM scale over which the transition relaxes). recombdist touches
# only the transition, so each swept value is EXACT. Truths are the SAME as the nnil figure:
#   SIM  = sim skim counts vs the simulated latent ancestry (BC2S3).
#   MolB = real skim counts vs nnil-on-molbreeding-hardcalls (the independent truth).
# Criterion = marker mismatch (argmin). Segment coordinates stay bp, so introgression
# sizes are in Mb exactly as for the other callers.
#
# Panels (2x2): A mismatch vs recombdist (log cM x), sim (blue) + MolB (orange), open
# circle = argmin; B sim ECDF at the two operating points vs latent ancestry; C real-skim
# ECDF vs MolB truth + ggtext KS; D QQ (sim y / molb x).
#
# Compute CACHED; NNIL_ZEALLB_RECOMPUTE=1 forces a rebuild.
#   Rscript scripts/fig_zeal_molbreeding_skim_lbimpute_calibration.R

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
RD_GRID <- c(1, 2, 5, 10, 20, 50, 100, 200) # recombdist in cM (transition relax scale; typical ~50)
R_MOLB <- 1.67e-3 # map estimate for the MolB truth caller (nnil on hardcalls)
TRUTH_CFG <- data.frame(nir = 0.7, germ = 1e-3, gert = 1e-4, p = 0.9, mr = 0)
KEEP <- c("name", "chr", "start_bp", "end_bp", "state")
BASE <- 20
ANN <- BASE * 0.8 / .pt
CACHE <- file.path(OUT, "fig_zeal_molbreeding_skim_lbimpute_calibration_cache.rds")

if (file.exists(CACHE) && Sys.getenv("NNIL_ZEALLB_RECOMPUTE") == "") {
  log_info("[lb] reusing cached compute (%s)", basename(CACHE))
  D <- readRDS(CACHE)
} else {
  seg_at <- function(seg, v) seg[abs(recombdist - v) < 1e-9]
  keyv <- function(v) sprintf("%g", v)
  cm_map <- fread(here::here("data/zeal/markers_snp50k_cm.tsv"))[, .(chr = as.integer(chr), pos = as.integer(pos), cm = as.numeric(cm))]

  # ---- SIM leg: lbimpute on sim skim COUNTS in cM (native map); truth = latent ancestry ----
  sim <- readRDS(file.path(OUT, "zeal_nil_bc2s3_full.rds"))
  N_SIM <- as.integer(Sys.getenv("NNIL_ZEALLB_NSIM", "500"))
  set.seed(1)
  sub <- sort(sample(seq_along(sim$names), min(N_SIM, length(sim$names))))
  sub_names <- sim$names[sub]
  M <- nrow(sim$grid)
  grid_cm <- data.table(chr = as.integer(sim$grid$chr), pos = as.integer(sim$grid$pos))
  grid_cm[cm_map, cm := i.cm, on = c("chr", "pos")] # attach the sim's own cM (order-preserving)
  if (anyNA(grid_cm$cm)) stop("sim markers missing cM in markers_snp50k_cm.tsv")
  sim_long <- data.table(
    name = rep(sub_names, each = M),
    chr = rep(grid_cm$chr, length(sub)),
    pos = rep(grid_cm$pos, length(sub)),
    cm = rep(grid_cm$cm, length(sub)),
    n_ref = as.integer(sim$n_ref[, sub]),
    n_alt = as.integer(sim$n_alt[, sub])
  )
  sim_truth <- as.data.table(sim$truth)[name %in% sub_names]
  sim_grid <- data.table(chr = as.integer(sim$grid$chr), pos = as.integer(sim$grid$pos))
  sim_traster <- rasterize_named(sim_truth, sim_grid)
  sim_tblocks <- .donor_blocks(sim_truth)
  threads <- min(parallel::detectCores() - 2L, 8L)
  log_info("[lb] SIM: %d/%d lines x %d markers; recombdist(cM) grid %s", length(sub), length(sim$names), M, paste(RD_GRID, collapse = ","))
  t0 <- Sys.time()
  sim_seg <- as.data.table(caller_sweep(sim_long, caller = "lbimpute", values = RD_GRID, unit = "cm", design = "BC2S3", threads = threads))
  log_info("[lb] SIM lbimpute sweep done in %.0fs", as.numeric(difftime(Sys.time(), t0, units = "secs")))
  sim_mm <- sapply(RD_GRID, function(v) 1 - marker_dsc(seg_at(sim_seg, v)[, ..KEEP], sim_truth, sim_grid, truth_raster = sim_traster)$accuracy)
  sim_dsc <- sapply(RD_GRID, function(v) donor_fragment_dsc(seg_at(sim_seg, v)[, ..KEEP], sim_truth, truth_blocks = sim_tblocks)$dsc)

  # ---- MolB truth (nnil on molbreeding hard calls); skim test: lbimpute on COUNTS in cM ----
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
  skim <- merge(skim, cm_map, by = c("chr", "pos"))
  if (anyNA(skim$cm)) stop("skim markers missing cM in markers_snp50k_cm.tsv")
  setorder(skim, name, chr, pos)
  skim_seg <- as.data.table(caller_sweep(skim[, .(name, chr, pos, cm, n_ref, n_alt)],
    caller = "lbimpute", values = RD_GRID, unit = "cm", design = "BC2S3", threads = 4L
  ))
  skim_mm <- sapply(RD_GRID, function(v) 1 - marker_dsc(seg_at(skim_seg, v)[, ..KEEP], molb_truth, molb_grid, truth_raster = molb_traster)$accuracy)
  skim_dsc <- sapply(RD_GRID, function(v) donor_fragment_dsc(seg_at(skim_seg, v)[, ..KEEP], molb_truth, truth_blocks = molb_tblocks)$dsc)

  D <- list(
    sweep = rbind(
      data.table(truth = "sim", recombdist = RD_GRID, mismatch = sim_mm, dsc = sim_dsc),
      data.table(truth = "molb", recombdist = RD_GRID, mismatch = skim_mm, dsc = skim_dsc)
    ),
    sim_truth_sz = donor_block_sizes(sim_truth),
    molb_truth_sz = donor_block_sizes(molb_truth),
    sim_sz_by_rd = setNames(lapply(RD_GRID, function(v) donor_block_sizes(seg_at(sim_seg, v)[, ..KEEP])), keyv(RD_GRID)),
    skim_sz_by_rd = setNames(lapply(RD_GRID, function(v) donor_block_sizes(seg_at(skim_seg, v)[, ..KEEP])), keyv(RD_GRID))
  )
  saveRDS(D, CACHE)
  fwrite(D$sweep, file.path(OUT, "fig_zeal_molbreeding_skim_lbimpute_calibration_sweep.csv"))
}
list2env(D, environment())

# ============================== figure ==============================
col_sim <- "#0072B2"
col_molb <- "#D55E00"
MB_BREAKS <- c(0.1, 1, 10, 100)
keyv <- function(v) sprintf("%g", v)

rd_sim <- sweep[truth == "sim"]$recombdist[which.min(sweep[truth == "sim"]$mismatch)]
rd_molb <- sweep[truth == "molb"]$recombdist[which.min(sweep[truth == "molb"]$mismatch)]
log_info("[lb] mismatch operating points: sim recombdist=%g cM | molb recombdist=%g cM", rd_sim, rd_molb)

sim_sz_sim <- sim_sz_by_rd[[keyv(rd_sim)]]
sim_sz_molb <- sim_sz_by_rd[[keyv(rd_molb)]]
skim_sz_sim <- skim_sz_by_rd[[keyv(rd_sim)]]
skim_sz_molb <- skim_sz_by_rd[[keyv(rd_molb)]]

L_sim <- sprintf("sim calibrated \nlbimpute (rd=%g cM)", rd_sim)
L_molb <- sprintf("molb calibrated \nlbimpute (rd=%g cM)", rd_molb)
L_st <- "sim latent ancestry\nBC2S3"
L_mt <- "molb calls"
pal <- c(setNames(c(col_sim, col_molb), c(L_sim, L_molb)), setNames(c("black", "black"), c(L_st, L_mt)))
lty <- c(setNames(c("solid", "solid"), c(L_sim, L_molb)), setNames(c("dotted", "dotted"), c(L_st, L_mt)))
allsz <- c(sim_truth_sz, sim_sz_sim, sim_sz_molb, molb_truth_sz, skim_sz_sim, skim_sz_molb)
xlim_mb <- c(0.1, max(allsz[is.finite(allsz) & allsz > 0]))

# A: mismatch vs recombdist (log cM x); open circle = mismatch-opt
optrow <- sweep[, .SD[which.min(mismatch)], by = truth]
p_A <- ggplot(sweep, aes(recombdist, mismatch, colour = truth)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.4) +
  geom_point(data = optrow, size = 3.6, shape = 21, fill = "white", stroke = 1.1) +
  scale_x_log10() +
  scale_colour_manual(
    values = c(sim = col_sim, molb = col_molb),
    labels = c(sim = "sim skim vs sim ancestry", molb = "real skim vs molb ancestry"), name = "calibration"
  ) +
  labs(x = "lbimpute recombdist (cM)", y = "marker mismatch rate", title = "lbimpute recombdist calibration") +
  theme_bw(base_size = BASE) +
  theme(
    aspect.ratio = 1, legend.position = c(0.02, 0.99), legend.justification = c(0, 1),
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
p_B <- mk_ecdf(ecdf_B, "Simulation:\nlbimpute calls vs latent ancestry", c(L_st, L_sim, L_molb))
p_C <- mk_ecdf(ecdf_C, "Real skim:\nlbimpute calls vs MolBreeding truth", c(L_mt, L_molb, L_sim))

ks_C <- suppressWarnings(ks.test(skim_sz_sim, skim_sz_molb))
ks_txt <- if (ks_C$p.value >= 0.01) sprintf("%.2f", ks_C$p.value) else sprintf("%.0e", ks_C$p.value)
log_info("[lb] panel C KS (sim-cal vs molb-cal lbimpute): D=%.3f p=%.3g", ks_C$statistic, ks_C$p.value)
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
    title = "QQ: sim- vs MolBreeding-\ncalibrated lbimpute (real skim)"
  ) +
  theme_bw(base_size = BASE) +
  theme(aspect.ratio = 1)

fig <- (p_A | p_B) / (p_C | p_D) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 25, face = "bold"),
    plot.tag.location = "plot", plot.tag.position = "topleft"
  )
ggsave(file.path(OUT, "fig_zeal_molbreeding_skim_lbimpute_calibration.png"), fig, width = 13, height = 13, dpi = 150)
log_info("[lb] wrote fig_zeal_molbreeding_skim_lbimpute_calibration.png")
