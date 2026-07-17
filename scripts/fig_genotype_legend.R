#!/usr/bin/env Rscript
# Standalone palette legend for the SNP50K genotype raster:
# REF gold, HET green, ALT purple, missing grey. Place beside the raster.

suppressPackageStartupMessages(library(ggplot2))

GOLD <- "#E0A81C"
GREEN <- "#2E9B57"
PURPLE <- "#6B3FA0"
GREY <- "#9AA0A6"
INK <- "#1A1A1A"

leg <- data.frame(
  x   = c(0, 1.4, 2.8, 4.2),
  lab = c("REF", "HET", "ALT", "MIS"),
  col = c(GOLD, GREEN, PURPLE, GREY)
)

gl <- ggplot(leg) +
  geom_tile(aes(x = x, y = 0),
    fill = leg$col, width = 0.4, height = 0.7,
    color = "grey30", linewidth = 0.5
  ) +
  geom_text(aes(x = x + 0.3, y = 0, label = lab), hjust = 0, size = 9, color = INK) +
  scale_x_continuous(limits = c(-0.35, 5.3)) +
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  theme_void()

out <- "nilhmm-paper/figures/genotype_raster_legend.png"
ggsave(out, gl, width = 6.6, height = 1.1, dpi = 150, bg = "white")
cat("wrote:", out, "\n")

## two-state legend: REF gold, ALT/HET purple
leg2 <- data.frame(x = c(0, 1.7), lab = c("REF", "ALT/HET"), col = c(GOLD, PURPLE))
gl2 <- ggplot(leg2) +
  geom_tile(aes(x = x, y = 0),
    fill = leg2$col, width = 0.4, height = 0.7,
    color = "grey30", linewidth = 0.5
  ) +
  geom_text(aes(x = x + 0.3, y = 0, label = lab), hjust = 0, size = 9, color = INK) +
  scale_x_continuous(limits = c(-0.35, 3.7)) +
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  theme_void()

out2 <- "nilhmm-paper/figures/genotype_raster_legend_2state.png"
ggsave(out2, gl2, width = 4.4, height = 1.1, dpi = 150, bg = "white")
cat("wrote:", out2, "\n")
