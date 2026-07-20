#!/usr/bin/env Rscript
# Calibration foil, step 6: the overlay figure.
#
# A  Calibration curves (donor-fragment Dice vs rrate) for the two truth sources:
#    the DENSE simcross truth (sharp interior optimum at rrate_sim*) and the SPARSE
#    24-line chip calls (no interior optimum -- monotone). Holland's chip-selected
#    avg_r and rrate_sim* are marked; the chip mismatch-admissible band is shaded.
#    Message: the simulation resolves the operating point that the chip leaves
#    broadly underdetermined, and lands on Holland's pick.
# B  Holland's own objective on the chip side (GBS-vs-chip calls mismatch): flat
#    across ~4 orders of magnitude of rrate -- WHY the chip cannot pin rrate down.
#
# Reads the step-4/5 sweep CSVs; writes nilhmm-paper/figures/nnil_calibration_foil.png.
#   Rscript scripts/nnil_foil/05_foil_figure.R

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})
root <- here::here()
FOIL <- file.path(root, "data/nnil_foil")
OUT <- file.path(root, "nilhmm-paper/figures/nnil_calibration_foil.png")

chip <- fread(file.path(FOIL, "chip_rrate_sweep.csv"))
sim <- fread(file.path(FOIL, "sim_rrate_sweep.csv"))
cj <- fromJSON(file.path(FOIL, "chip_calib.json"))

map_r <- cj$map_r # map-defined per-marker recombination fraction (native v5 map)
rrate_sim <- sim$rrate[which.max(sim$donor_frag_dice)]
band <- cj$mismatch_plateau # chip-admissible rrate band (flat mismatch)

col_sim <- "#0072B2" # Okabe-Ito blue
col_chip <- "#D55E00" # Okabe-Ito vermillion
col_holl <- "grey30"

base <- theme_classic(base_size = 9) +
  theme(
    text = element_text(family = "sans"),
    plot.tag = element_text(face = "bold"),
    legend.position = "top", legend.title = element_blank(),
    legend.key.height = unit(0.4, "lines"),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3)
  )
lx <- scale_x_log10(
  breaks = 10^(-6:-1),
  labels = c("1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1")
)

# ---- Panel A: calibration curves --------------------------------------------
cur <- rbind(
  data.table(rrate = sim$rrate, dice = sim$donor_frag_dice, src = "Simulation (dense truth)"),
  data.table(rrate = chip$rrate, dice = chip$donor_frag_dice, src = "Chip (24-line truth)")
)
pA <- ggplot(cur, aes(rrate, dice, colour = src)) +
  annotate("rect",
    xmin = band[1], xmax = band[2], ymin = -Inf, ymax = Inf,
    fill = "grey85", alpha = 0.5
  ) +
  geom_vline(xintercept = map_r, linetype = "dashed", colour = col_holl, linewidth = 0.4) +
  geom_vline(xintercept = rrate_sim, linetype = "dotted", colour = col_sim, linewidth = 0.5) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1) +
  scale_colour_manual(values = c("Simulation (dense truth)" = col_sim, "Chip (24-line truth)" = col_chip)) +
  lx +
  annotate("text",
    x = map_r, y = 0.30, label = "map~italic(r)", parse = TRUE,
    angle = 90, vjust = -0.4, hjust = 0, size = 2.6, colour = col_holl
  ) +
  annotate("text",
    x = rrate_sim, y = 0.98, label = "italic(r)[sim]^'*'", parse = TRUE,
    angle = 90, vjust = 1.3, hjust = 1, size = 2.6, colour = col_sim
  ) +
  labs(
    x = expression("intermarker recombination fraction, " * italic(r)),
    y = "Donor-fragment DSC", tag = "A"
  ) +
  base
ann_band <- sprintf("chip-admissible band\n(flat mismatch, %.0e to %.0e)", band[1], band[2])
pA <- pA + annotate("text",
  x = sqrt(band[1] * band[2]), y = 0.12,
  label = "chip-admissible band", size = 2.5, colour = "grey40"
)

# ---- Panel B: chip mismatch is flat -----------------------------------------
pB <- ggplot(chip, aes(rrate, holland_mismatch)) +
  annotate("rect",
    xmin = band[1], xmax = band[2], ymin = -Inf, ymax = Inf,
    fill = "grey85", alpha = 0.5
  ) +
  geom_vline(xintercept = map_r, linetype = "dashed", colour = col_holl, linewidth = 0.4) +
  geom_line(linewidth = 0.6, colour = col_chip) +
  geom_point(size = 1, colour = col_chip) +
  lx +
  annotate("text",
    x = map_r, y = mean(range(chip$holland_mismatch)), label = "map~italic(r)",
    parse = TRUE, angle = 90, vjust = -0.4, hjust = 0.5, size = 2.6, colour = col_holl
  ) +
  labs(
    x = expression("intermarker recombination fraction, " * italic(r)),
    y = "GBS-vs-chip calls mismatch", tag = "B"
  ) +
  base +
  theme(legend.position = "none")

fig <- pA + pB + plot_layout(widths = c(1.25, 1))
ggsave(OUT, fig, width = 180, height = 78, units = "mm", dpi = 300)
cat(sprintf(
  "wrote %s\n  rrate_sim* = %.3e | map r = %.3e | chip band [%.1e, %.1e]\n",
  OUT, rrate_sim, map_r, band[1], band[2]
))
