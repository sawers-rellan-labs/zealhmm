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
variant <- commandArgs(TRUE)[2]
if (is.na(variant)) variant <- "ols" # "ols" or "mlm"
sfx <- if (variant == "mlm") "_mlm" else ""
mlab <- if (variant == "mlm") "MLM (Family+K)" else "OLS"
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

sw <- as.data.table(fread(sprintf("results/sim/teonam/stam_gwas_%s_118k%s_sweep.csv", caller, sfx)))
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
ytop <- ceiling(max(sw$logP)) # automatic y-limit (fit to the data)
sw[, dord := ifelse(is.finite(coverage), -coverage, Inf)]

# LOD-5 peak loci (above-threshold markers clumped within 1 Mb -> top marker per clump)
# to annotate: infinity = black up-triangles just above the axis; 0.5x = blue
# down-triangles just below the top (mirror image).
peak_loci <- function(cvsub, gap = 1e6) {
  s <- cvsub[logP > LOD]
  if (!nrow(s)) {
    return(s[0])
  }
  setorder(s, CHR, BP)
  s[, pk := cumsum(CHR != shift(CHR, fill = -1L) | BP - shift(BP, fill = -Inf) > gap)]
  s[, .SD[which.max(logP)], by = pk]
}
pk_inf <- peak_loci(sw[!is.finite(coverage)])
pk_05 <- peak_loci(sw[coverage == 0.5])

p <- ggplot(sw[order(dord, CHR, BP)], aes(x = BPn, y = logP, color = cov_lab, group = interaction(CHR, coverage))) +
  geom_line(linewidth = 0.4) +
  geom_hline(yintercept = LOD, linetype = "dotted", linewidth = 0.7, color = "grey30") +
  geom_point(
    data = pk_inf, aes(x = BPn, y = -0.045 * ytop), inherit.aes = FALSE,
    shape = 24, fill = "black", colour = "black", size = 2.8
  ) + # infinity LOD-5 loci: black up-triangles BELOW the x-axis (shape 24 mirrors 25)
  geom_point(
    data = pk_05, aes(x = BPn, y = 0.98 * ytop), inherit.aes = FALSE,
    shape = 25, fill = pal[["0.5"]], colour = pal[["0.5"]], size = 2.8
  ) + # 0.5x LOD-5 loci: down-triangles in the 0.5x legend colour, at the top
  scale_color_manual(values = pal, name = "coverage (×)", drop = FALSE) +
  scale_x_continuous(breaks = axis_df$mid, labels = axis_df$CHR, expand = c(0.01, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0))) +
  coord_cartesian(ylim = c(0, ytop), clip = "off") + # allow the below-axis triangles to render in the margin
  labs(
    x = "chromosome", y = expression(-log[10](italic(P))),
    title = sprintf("STAM — %s ancestry, %s GWAS −log10P vs coverage (authentic 118K truth)", toupper(caller), mlab)
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.margin = margin(t = 5, r = 6, b = 16, l = 6), # room below the axis for the triangles
    plot.title = element_markdown(hjust = 0, size = 14, face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
    legend.position = c(0.995, 0.97), legend.justification = c(1, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.7), colour = NA),
    legend.key.size = unit(0.85, "lines"),
    legend.title = element_text(size = 12), legend.text = element_text(size = 11)
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.2)))

out <- sprintf("results/sim/teonam/stam_gwas_%s_118k%s_sweep_manhattan.png", caller, sfx)
ggsave(out, p, width = 9, height = 4.6, dpi = 200, bg = "white")
cat("wrote", out, "| ytop =", ytop, "\n")
