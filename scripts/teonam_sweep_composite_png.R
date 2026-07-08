#!/usr/bin/env Rscript
# Four-panel coverage-degradation composite (control / nNIL / RTIGER / LB-Impute) for
# one GWAS variant. Same plot_sweep_line as the notebook.
#   Rscript scripts/teonam_sweep_composite_png.R mlm   # MLM (Family+K) sweeps
#   Rscript scripts/teonam_sweep_composite_png.R ols   # OLS sweeps
suppressMessages({
  library(data.table)
  library(ggplot2)
  library(ggtext)
  library(scales)
  library(cowplot)
})
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)
variant <- commandArgs(TRUE)[1]
if (is.na(variant)) variant <- "mlm"
sfx <- if (variant == "mlm") "_mlm" else ""
LOD <- 5

get_transformer <- function(m) {
  cmat <- do.call(rbind, lapply(split(m[order(m$CHR, m$BP), ], m$CHR[order(m$CHR, m$BP)]), function(d) {
    data.frame(CHR = d$CHR[1], min_bp = min(d$BP), width = max(d$BP) - min(d$BP), medgap = if (nrow(d) > 1) median(diff(sort(d$BP))) else 1)
  }))
  cmat <- cmat[order(cmat$CHR), ]
  maxgap <- max(cmat$medgap, na.rm = TRUE)
  numc <- nrow(cmat)
  cmat$base <- 0
  cmat$midp <- 0
  cmat$midp[1] <- cmat$width[1] / 2
  for (i in 2:numc) {
    cmat$base[i] <- cmat$base[i - 1] + cmat$width[i - 1] + maxgap
    cmat$midp[i] <- cmat$base[i] + cmat$width[i] / 2
  }
  fac <- numc / cmat$midp[numc]
  cmat$basef <- fac * cmat$base
  function(chr, bp) fac * (bp - cmat$min_bp[match(chr, cmat$CHR)]) + cmat$basef[match(chr, cmat$CHR)]
}

plot_sweep_line <- function(csv, title, legend = TRUE, ytop = 20) {
  sw <- as.data.table(fread(csv))[is.finite(P) & P > 0]
  sw[, logP := -log10(P)]
  tr <- get_transformer(as.data.frame(sw[coverage == sw$coverage[1], .(CHR, BP)]))
  sw[, BPn := tr(CHR, BP)]
  lev <- c("∞", "20", "10", "5", "1", "0.5", "0.2", "0.1")
  sw[, cov_lab := factor(ifelse(is.finite(coverage), formatC(coverage), "∞"), levels = lev)]
  pal <- setNames(grDevices::hcl.colors(7, "viridis"), c("0.1", "0.2", "0.5", "1", "5", "10", "20"))
  pal["∞"] <- "black"
  axis_df <- sw[, .(mid = (min(BPn) + max(BPn)) / 2), by = CHR][order(CHR)]
  sw[, dord := ifelse(is.finite(coverage), -coverage, Inf)]
  ggplot(sw[order(dord, CHR, BP)], aes(BPn, logP, color = cov_lab, group = interaction(CHR, coverage))) +
    geom_line(linewidth = 0.4) +
    geom_hline(yintercept = LOD, linetype = "dotted", linewidth = 0.7, color = "grey30") +
    scale_color_manual(values = pal, name = "coverage (×)", drop = FALSE) +
    scale_x_continuous(breaks = axis_df$mid, labels = axis_df$CHR, expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, ytop), expand = expansion(mult = c(0.01, 0))) +
    labs(x = "chromosome", y = expression(-log[10](italic(P))), title = title) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_markdown(hjust = 0, size = 14, face = "bold"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
      legend.position = if (legend) c(0.995, 0.97) else "none", legend.justification = c(1, 1),
      legend.background = element_rect(fill = scales::alpha("white", 0.7), colour = NA),
      legend.key.size = unit(0.85, "lines"), legend.title = element_text(size = 11), legend.text = element_text(size = 10)
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 1.2)))
}

f <- function(c) file.path(ROOT, sprintf("results/sim/teonam/stam_gwas_%s_118k%s_sweep.csv", c, sfx))
lab <- if (variant == "mlm") "MLM (Family+K)" else "OLS"
pA <- plot_sweep_line(f("control"), sprintf("STAM %s — interpolation control (GL+HWE)", lab), legend = TRUE)
pB <- plot_sweep_line(f("nnil"), sprintf("STAM %s — nNIL ancestry", lab), legend = FALSE)
pC <- plot_sweep_line(f("rtiger"), sprintf("STAM %s — RTIGER ancestry", lab), legend = FALSE)
pD <- plot_sweep_line(f("lbimpute"), sprintf("STAM %s — LB-Impute ancestry", lab), legend = FALSE)
comp <- cowplot::plot_grid(pA, pB, pC, pD, ncol = 1, align = "v", axis = "lr", labels = c("A", "B", "C", "D"), label_size = 20)
out <- file.path(ROOT, sprintf("results/sim/teonam/stam_sweep_composite%s_118k.png", sfx))
ggsave(out, comp, width = 9, height = 16, dpi = 200, bg = "white")
cat("wrote", out, "\n")
