#!/usr/bin/env Rscript
# Does googa (real BrB, rrate=map) reproduce the fragment-size distribution of BOTH independent
# references: the BC2S3 SIM latent target AND the MolBreeding data (nnil-on-molb-hardcalls,
# different samples, distributional only)? Reuses the cached googa nir segments and both cached
# size distributions -- no decode.
# CAVEAT: molb sizes are on the wsfilt grid (~9k markers), sim/googa on the pangene grid, so
# cross-reference size differences partly reflect grid resolution, not just caller quality.
#   Rscript scripts/fig_zeal_googa_vs_molb.R

suppressMessages({
  library(data.table)
  library(here)
  library(ggplot2)
  library(patchwork)
})
devtools::load_all(path.expand("~/repos/nilhmm"), quiet = TRUE)
source(here::here("R/metrics.R"))
OUT <- here::here("results/sim/zeal_nil")
BASE <- 20
KEEP <- c("name", "chr", "start_bp", "end_bp", "state")

sim_sz <- readRDS(file.path(OUT, "atlas_target_sz.rds"))
molb_sz <- readRDS(file.path(OUT, "fig_zeal_molbreeding_skim_bbnil_calibration_cache.rds"))$molb_truth_sz
SEG_DIR <- file.path(OUT, "googa_nir_seg_cache")
files <- list.files(SEG_DIR, pattern = "^nir_.*\\.rds$", full.names = TRUE)
nirs <- as.numeric(sub("nir_(.*)\\.rds", "\\1", basename(files)))
o <- order(nirs)
files <- files[o]
nirs <- nirs[o]

sweep <- rbindlist(lapply(seq_along(files), function(i) {
  sz <- donor_block_sizes(as.data.table(readRDS(files[i]))[, ..KEEP])
  data.table(
    nir = nirs[i], ks_sim = fragment_size_ks(sz, sim_sz),
    ks_molb = fragment_size_ks(sz, molb_sz), median_mb = stats::median(sz)
  )
}))
cat(sprintf(
  "references: sim median=%.2f Mb (n=%d) | molb median=%.2f Mb (n=%d)\n",
  median(sim_sz), length(sim_sz), median(molb_sz), length(molb_sz)
))
print(sweep, digits = 3)
nir_sim <- sweep$nir[which.min(sweep$ks_sim)]
nir_molb <- sweep$nir[which.min(sweep$ks_molb)]
cat(sprintf(
  "KS-opt vs SIM: nir=%.2f (D=%.3f) | KS-opt vs MOLB: nir=%.2f (D=%.3f)\n",
  nir_sim, min(sweep$ks_sim), nir_molb, min(sweep$ks_molb)
))
fwrite(sweep, file.path(OUT, "fig_zeal_googa_vs_molb_sweep.csv"))

# ---- figure: KS-vs-nir for both refs (A) + ECDF overlay (B) ------------------
long <- melt(sweep,
  id.vars = "nir", measure.vars = c("ks_sim", "ks_molb"),
  variable.name = "ref", value.name = "ks"
)
long[, ref := factor(ref, c("ks_sim", "ks_molb"), c("vs BC2S3 sim", "vs MolBreeding"))]
p_A <- ggplot(long, aes(nir, ks, colour = ref)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c("vs BC2S3 sim" = "#0072B2", "vs MolBreeding" = "#CC79A7"), name = NULL) +
  expand_limits(y = 0) +
  labs(
    x = "googa nir (rrate pinned to map)", y = "introgression-size KS distance (D)",
    title = "googa (real BrB) fragment sizes\nvs sim and vs MolBreeding"
  ) +
  theme_bw(base_size = BASE) +
  theme(
    aspect.ratio = 1, legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
    legend.background = element_rect(fill = "transparent", colour = NA)
  )

sz_molbopt <- donor_block_sizes(as.data.table(readRDS(file.path(SEG_DIR, sprintf("nir_%.3f.rds", nir_molb))))[, ..KEEP])
ecdf_dt <- rbindlist(list(
  data.table(size_mb = molb_sz, series = "MolBreeding (nnil, wsfilt grid)"),
  data.table(size_mb = sim_sz, series = "BC2S3 sim latent (pangene grid)"),
  data.table(size_mb = sz_molbopt, series = sprintf("googa BrB (nir=%.2f)", nir_molb))
))
p_B <- ggplot(ecdf_dt, aes(size_mb, colour = series)) +
  stat_ecdf(linewidth = 1) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100)) +
  scale_colour_manual(values = setNames(
    c("#CC79A7", "black", "#0072B2"),
    c("MolBreeding (nnil, wsfilt grid)", "BC2S3 sim latent (pangene grid)", sprintf("googa BrB (nir=%.2f)", nir_molb))
  ), name = NULL) +
  labs(x = "introgression size (Mb)", y = "ECDF", title = "Fragment-size distributions") +
  theme_bw(base_size = BASE) +
  theme(
    aspect.ratio = 1, legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
    legend.background = element_rect(fill = "transparent", colour = NA)
  )

fig <- (p_A | p_B) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 22, face = "bold"))
ggsave(file.path(OUT, "fig_zeal_googa_vs_molb.png"), fig, width = 12, height = 6, dpi = 150)
cat("wrote fig_zeal_googa_vs_molb.png\n")
