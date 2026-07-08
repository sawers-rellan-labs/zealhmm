#!/usr/bin/env Rscript
# Single-coverage Manhattan (point plot) for one caller/variant/lambda from a sweep CSV.
#   Rscript scripts/teonam_single_manhattan_png.R rtiger mlm 1
suppressMessages({
  library(data.table)
  library(ggplot2)
  library(ggtext)
})
ROOT <- "/Users/fvrodriguez/repos/zealhmm"
setwd(ROOT)
a <- commandArgs(TRUE)
caller <- if (length(a) >= 1) a[1] else "rtiger"
variant <- if (length(a) >= 2) a[2] else "mlm"
lam <- if (length(a) >= 3) as.numeric(a[3]) else 1
sfx <- if (variant == "mlm") "_mlm" else ""
mlab <- if (variant == "mlm") "MLM (Family+K)" else "OLS"
LOD <- 5

get_transformer <- function(m) {
  cmat <- do.call(rbind, lapply(split(m[order(m$CHR, m$BP), ], m$CHR[order(m$CHR, m$BP)]), function(d) {
    data.frame(CHR = d$CHR[1], min_bp = min(d$BP), width = max(d$BP) - min(d$BP))
  }))
  cmat <- cmat[order(cmat$CHR), ]
  numc <- nrow(cmat)
  gap <- 0.02 * sum(cmat$width) / numc
  cmat$base <- 0
  for (i in 2:numc) cmat$base[i] <- cmat$base[i - 1] + cmat$width[i - 1] + gap
  function(chr, bp) cmat$base[match(chr, cmat$CHR)] + (bp - cmat$min_bp[match(chr, cmat$CHR)])
}

sw <- fread(sprintf("results/sim/teonam/stam_gwas_%s_118k%s_sweep.csv", caller, sfx))
s <- sw[coverage == lam & is.finite(P) & P > 0]
s[, logP := -log10(P)]
setorder(s, CHR, BP)
tr <- get_transformer(as.data.frame(s[, .(CHR, BP)]))
s[, x := tr(CHR, BP)]
s[, oddeven := ifelse(CHR %% 2 == 0, "e", "o")]
axis_df <- s[, .(mid = (min(x) + max(x)) / 2), by = CHR][order(CHR)]
covlab <- if (is.infinite(lam)) "∞" else as.character(lam)
ytop <- max(20, ceiling(max(s$logP)))

p <- ggplot(s, aes(x, logP, color = oddeven)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = LOD, linetype = "dotted", color = "grey30") +
  scale_color_manual(values = c(o = "black", e = "grey65"), guide = "none") +
  scale_x_continuous(breaks = axis_df$mid, labels = axis_df$CHR, expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, ytop), expand = expansion(mult = c(0.01, 0))) +
  labs(
    x = "chromosome", y = expression(-log[10](italic(P))),
    title = sprintf("STAM — %s ancestry, %s GWAS, coverage = %s× (authentic 118K truth)", toupper(caller), mlab, covlab)
  ) +
  theme_classic(base_size = 15) +
  theme(plot.title = element_markdown(size = 13, face = "bold"), panel.grid.major.y = element_line(color = "grey92", linewidth = 0.2))

out <- sprintf("results/sim/teonam/stam_gwas_%s_118k%s_lambda%s_manhattan.png", caller, sfx, lam)
ggsave(out, p, width = 9, height = 4.4, dpi = 200, bg = "white")
cat("wrote", out, "| ytop =", ytop, "| markers =", nrow(s), "\n")
