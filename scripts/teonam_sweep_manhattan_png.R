#!/usr/bin/env Rscript
# Standalone line-Manhattan PNG for one 118K caller sweep (same plot_sweep_line as
# analysis/teonam-qtl-recovery-118k.qmd). Reads the combined sweep CSV, draws
# -log10P vs genome position, one line per coverage (viridis; lambda=Inf in black).
# Run: Rscript scripts/teonam_sweep_manhattan_png.R rtiger
suppressMessages({
  library(data.table)
  library(ggplot2)
  library(ggtext)
  library(scales)
})
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)
caller <- commandArgs(TRUE)[1]
if (is.na(caller)) caller <- "rtiger"
LOD <- 5

get_transformer <- function(m) {
  cmat <- do.call(rbind, lapply(
    split(m[order(m$CHR, m$BP), ], m$CHR[order(m$CHR, m$BP)]),
    function(d) {
      data.frame(
        CHR = d$CHR[1], min_bp = min(d$BP), width = max(d$BP) - min(d$BP),
        medgap = if (nrow(d) > 1) median(diff(sort(d$BP))) else 1
      )
    }
  ))
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
  function(chr, bp) {
    i <- match(chr, cmat$CHR)
    fac * (bp - cmat$min_bp[i]) + cmat$basef[i]
  }
}

sw <- as.data.table(fread(sprintf("results/sim/teonam/stam_gwas_%s_118k_sweep.csv", caller)))
sw <- sw[is.finite(P) & P > 0]
sw[, logP := -log10(P)]
tr <- get_transformer(as.data.frame(sw[coverage == sw$coverage[1], .(CHR, BP)]))
sw[, BPn := tr(CHR, BP)]
cov_levels <- c("∞", "20", "10", "5", "1", "0.5", "0.2", "0.1")
sw[, cov_lab := factor(ifelse(is.finite(coverage), formatC(coverage), "∞"), levels = cov_levels)]
cov_asc <- c("0.1", "0.2", "0.5", "1", "5", "10", "20")
pal <- setNames(grDevices::hcl.colors(7, "viridis"), cov_asc)
pal["∞"] <- "black"
axis_df <- sw[, .(mid = (min(BPn) + max(BPn)) / 2), by = CHR][order(CHR)]
ytop <- max(20, ceiling(max(sw$logP))) # adaptive so strong peaks are not clipped
sw[, dord := ifelse(is.finite(coverage), -coverage, Inf)]

p <- ggplot(sw[order(dord, CHR, BP)], aes(x = BPn, y = logP, color = cov_lab, group = interaction(CHR, coverage))) +
  geom_line(linewidth = 0.4) +
  geom_hline(yintercept = LOD, linetype = "dotted", linewidth = 0.7, color = "grey30") +
  scale_color_manual(values = pal, name = "coverage (×)", drop = FALSE) +
  scale_x_continuous(breaks = axis_df$mid, labels = axis_df$CHR, expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, ytop), expand = expansion(mult = c(0.01, 0))) +
  labs(
    x = "chromosome", y = expression(-log[10](italic(P))),
    title = sprintf("STAM — %s ancestry, GWAS −log10P vs coverage (authentic 118K truth)", toupper(caller))
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_markdown(hjust = 0, size = 14, face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
    legend.position = c(0.995, 0.97), legend.justification = c(1, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.7), colour = NA),
    legend.key.size = unit(0.85, "lines"),
    legend.title = element_text(size = 12), legend.text = element_text(size = 11)
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.2)))

out <- sprintf("results/sim/teonam/stam_gwas_%s_118k_sweep_manhattan.png", caller)
ggsave(out, p, width = 9, height = 4.6, dpi = 200, bg = "white")
cat("wrote", out, "| ytop =", ytop, "\n")
