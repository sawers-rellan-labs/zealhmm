#!/usr/bin/env Rscript
# Figure: genotype-call counts by read depth, one histogram per state
# (REF / HET / ALT), each with its OWN linear y-axis (facet scales = free_y).
# Shows that HET calls pile up at depth 1 (single-read -> prior forces HET) while
# ALT-hom stays rare at every depth.
#
#   depth   REF          HET       ALT
#   1       15,868,677   448,805   355
#   2       4,057,642    154,596   138
#   3-5     1,184,345    44,361    343
#   6-10    17,559       727       163
#   >10     173          49        2
# depth 0 dropped: those were "calls" with zero read support (REF 40,524, HET 1),
# which are not real calls; the 51.5M truly-missing sites are not plotted either.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtext)
  library(scales)
})

here <- tryCatch(dirname(sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)
)), error = function(e) ".")
if (length(here) == 0 || here == "") here <- "scripts"
FIGDIR <- normalizePath(file.path(here, "..", "nilhmm-paper", "figures"), mustWork = FALSE)

GOLD <- "#B8860B"
GREEN <- "#2E9B57"
PURPLE <- "#6B3FA0"
INK <- "#1A1A1A"
# (REF darkened from the paint gold #E0A81C so the strip label reads on white)

lvl <- c("1", "2", "3-5", "6-10", ">10")
df <- data.frame(
  depth = factor(rep(lvl, 3), levels = lvl),
  geno = factor(rep(c("REF", "HET", "ALT"), each = 5), levels = c("REF", "HET", "ALT")),
  count = c(
    15868677, 4057642, 1184345, 17559, 173, # REF
    448805, 154596, 44361, 727, 49, # HET
    355, 138, 343, 163, 2
  )
) # ALT
dfp <- df[df$count > 0, ] # log scale: drop zeros

strip <- as_labeller(c(
  REF = sprintf("<span style='color:%s'>**REF-hom**</span>", GOLD),
  HET = sprintf("<span style='color:%s'>**HET**</span>", GREEN),
  ALT = sprintf("<span style='color:%s'>**ALT-hom**</span>", PURPLE)
))

p <- ggplot(dfp, aes(depth, count, fill = geno)) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.2) +
  geom_text(aes(label = comma(count)), vjust = -0.45, size = 2.9, color = INK) +
  facet_wrap(~geno, scales = "free_y", labeller = strip) +
  scale_fill_manual(values = c(REF = GOLD, HET = GREEN, ALT = PURPLE), guide = "none") +
  scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "read depth", y = "calls",
    title = "Genotype calls by read depth",
    subtitle = "HET calls pile up at a single read (depth 1); ALT-hom stays rare at every depth."
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_markdown(size = 14),
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text = element_text(color = INK),
    axis.text.x = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 15, color = INK),
    plot.subtitle = element_text(size = 10.5, color = "grey30"),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 12, 8, 10)
  )

out <- file.path(FIGDIR, "calls_by_depth_histograms.png")
ggsave(out, p, width = 10, height = 4.2, dpi = 200, bg = "white")
cat("wrote:", out, "\n")
