#!/usr/bin/env Rscript
# Full 4-panel googa BrB calibration figure (rrate pinned to the pangene-grid map rate; nir the
# knob). TWO calibrations of the same googa-on-real-BrB run:
#   sim-calibrated  = nir minimizing KS to the BC2S3 sim-latent fragment sizes (pangene grid);
#   molb-calibrated = nir minimizing KS to the MolBreeding fragment sizes (nnil, wsfilt grid).
# A: KS vs nir for both references (both optima marked).
# B: BC2S3 sim latent (dotted) vs the two calibrations.
# C: MolBreeding (dotted) vs the two calibrations.
# D: QQ of the two calibrations against each other.
# Reuses cached googa nir segments + both reference size distributions -- no decode.
# CAVEAT: molb sizes are on the coarser wsfilt grid (~9k) vs sim/googa on the pangene grid (~23k);
# part of the molb-vs-sim size gap is grid resolution + nnil smoothing, not biology.
#   Rscript scripts/fig_zeal_googa_calibration.R

suppressMessages({
  library(data.table)
  library(here)
  library(ggplot2)
  library(patchwork)
  library(ggtext)
})
devtools::load_all(path.expand("~/repos/nilhmm"), quiet = TRUE)
source(here::here("R/metrics.R"))
OUT <- here::here("results/sim/zeal_nil")
BASE <- 20
KEEP <- c("name", "chr", "start_bp", "end_bp", "state")
SEG_DIR <- file.path(OUT, "googa_nir_seg_cache")
MB_BREAKS <- c(0.1, 1, 10, 100)

sim_sz <- readRDS(file.path(OUT, "atlas_target_sz.rds"))
molb_sz <- readRDS(file.path(OUT, "fig_zeal_molbreeding_skim_bbnil_calibration_cache.rds"))$molb_truth_sz
sweep <- unique(fread(file.path(OUT, "fig_zeal_googa_vs_molb_sweep.csv")), by = "nir")[order(nir)]
nir_sim <- sweep$nir[which.min(sweep$ks_sim)]
nir_molb <- sweep$nir[which.min(sweep$ks_molb)]
gsz <- function(v) donor_block_sizes(as.data.table(readRDS(file.path(SEG_DIR, sprintf("nir_%.3f.rds", v))))[, ..KEEP])
sz_simcal <- gsz(nir_sim)
sz_molbcal <- gsz(nir_molb)

L_sc <- sprintf("googa BrB sim-cal (nir=%.2f)", nir_sim)
L_mc <- sprintf("googa BrB molb-cal (nir=%.2f)", nir_molb)
L_sr <- "BC2S3 sim latent (pangene grid)"
L_mr <- "MolBreeding nnil (wsfilt grid)"
col_sc <- "#0072B2"
col_mc <- "#D55E00"
allsz <- c(sim_sz, molb_sz, sz_simcal, sz_molbcal)
xlim_mb <- c(0.05, max(allsz[is.finite(allsz) & allsz > 0]))

# ---- A: calibration curves --------------------------------------------------
long <- melt(sweep, id.vars = "nir", measure.vars = c("ks_sim", "ks_molb"), variable.name = "ref", value.name = "ks")
long[, ref := factor(ref, c("ks_sim", "ks_molb"), c("vs BC2S3 sim", "vs MolBreeding"))]
opt <- rbind(
  data.table(nir = nir_sim, ks = min(sweep$ks_sim), ref = "vs BC2S3 sim"),
  data.table(nir = nir_molb, ks = min(sweep$ks_molb), ref = "vs MolBreeding")
)
opt[, ref := factor(ref, levels = levels(long$ref))]
p_A <- ggplot(long, aes(nir, ks, colour = ref)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.6) +
  geom_point(data = opt, size = 3.6, shape = 21, fill = "white", stroke = 1.1) +
  scale_colour_manual(values = c("vs BC2S3 sim" = col_sc, "vs MolBreeding" = col_mc), name = "calibration target") +
  expand_limits(y = 0) +
  labs(
    x = "googa nir (rrate pinned to map)", y = "introgression-size KS (D)",
    title = "googa nir calibration:\nsim vs MolBreeding targets"
  ) +
  theme_bw(base_size = BASE) +
  theme(
    aspect.ratio = 1, legend.position = c(0.5, 0.98), legend.justification = c(0.5, 1),
    legend.background = element_rect(fill = "transparent", colour = NA)
  )

mk_ecdf <- function(ref_sz, ref_lab, ttl) {
  d <- rbindlist(list(
    data.table(size_mb = ref_sz, series = ref_lab, kind = "ref"),
    data.table(size_mb = sz_simcal, series = L_sc, kind = "cal"),
    data.table(size_mb = sz_molbcal, series = L_mc, kind = "cal")
  ))
  d[, series := factor(series, levels = c(ref_lab, L_sc, L_mc))]
  ggplot(d, aes(size_mb, colour = series, linetype = series)) +
    stat_ecdf(linewidth = 1) +
    scale_x_log10(breaks = MB_BREAKS, limits = xlim_mb, oob = scales::oob_keep) +
    scale_colour_manual(values = setNames(c("black", col_sc, col_mc), c(ref_lab, L_sc, L_mc)), name = NULL) +
    scale_linetype_manual(values = setNames(c("dotted", "solid", "solid"), c(ref_lab, L_sc, L_mc)), name = NULL) +
    labs(x = "introgression size (Mb)", y = "ECDF", title = ttl) +
    theme_bw(base_size = BASE) +
    theme(
      aspect.ratio = 1, legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
      legend.text = element_text(size = 10),
      legend.background = element_rect(fill = "transparent", colour = NA)
    )
}
p_B <- mk_ecdf(sim_sz, L_sr, "Simulation target\nvs the two calibrations")
p_C <- mk_ecdf(molb_sz, L_mr, "MolBreeding target\nvs the two calibrations")

# ---- D: QQ of the two calibrations against each other -----------------------
qs <- ppoints(200)
qq <- data.table(simcal = quantile(sz_simcal, qs), molbcal = quantile(sz_molbcal, qs))
ks_cal <- suppressWarnings(ks.test(sz_simcal, sz_molbcal))
p_D <- ggplot(qq, aes(molbcal, simcal)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50") +
  geom_point(size = 1.3, colour = col_sc) +
  scale_x_log10(limits = xlim_mb, breaks = MB_BREAKS, oob = scales::oob_keep) +
  scale_y_log10(limits = xlim_mb, breaks = MB_BREAKS, oob = scales::oob_keep) +
  labs(
    x = sprintf("size (Mb) molb-cal (nir=%.2f)", nir_molb),
    y = sprintf("size (Mb) sim-cal (nir=%.2f)", nir_sim),
    title = sprintf("QQ: sim- vs molb-\ncalibrated googa (D=%.2f)", ks_cal$statistic)
  ) +
  theme_bw(base_size = BASE) +
  theme(aspect.ratio = 1)

fig <- (p_A | p_B) / (p_C | p_D) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 25, face = "bold"), plot.tag.location = "plot", plot.tag.position = "topleft")
ggsave(file.path(OUT, "fig_zeal_googa_calibration.png"), fig, width = 13, height = 13, dpi = 150)
cat(sprintf(
  "sim-cal nir=%.2f | molb-cal nir=%.2f | QQ D=%.3f\nwrote fig_zeal_googa_calibration.png\n",
  nir_sim, nir_molb, ks_cal$statistic
))
