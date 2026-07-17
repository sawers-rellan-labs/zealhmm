#!/usr/bin/env Rscript
# Figure: genotype-call composition of the skim SNP50K matrix, as a progressive
# "zoom" stacked bar (technique after stefan, https://stackoverflow.com/a/76094965,
# CC BY-SA 4.0): each bar magnifies a sub-region of the one to its left, connected
# by grey trapezoids.
#
#   All calls  ->  Called (drop missing)  ->  Non-REF called (drop REF)
#
# Why three steps: ALT-hom is 0.00% of all AND 0.00% of called, so only the last
# zoom (into the non-REF calls) makes it visible -- exposing the het-excess:
# 648,539 HET vs 1,001 ALT-hom (~648x), the depth-1 calling artifact.
#
#   category        count        % of all   % of called
#   REF-hom 0/0     21,128,396   28.81      97.02
#   HET             648,538       0.88       2.98
#   ALT-hom          1,001        0.00       0.00
#   missing ./.     51,563,345   70.31      —
# (depth-0 "calls" with zero read support -- REF 40,524, HET 1 -- filtered out.)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtext)
})

here <- tryCatch(dirname(sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)
)), error = function(e) ".")
if (length(here) == 0 || here == "") here <- "scripts"
FIGDIR <- normalizePath(file.path(here, "..", "nilhmm-paper", "figures"), mustWork = FALSE)

GOLD <- "#E0A81C"
GREEN <- "#2E9B57"
PURPLE <- "#6B3FA0"
GREY <- "#9AA0A6"
INK <- "#1A1A1A"
RED <- "#C0392B"

## ---- counts ---------------------------------------------------------------
n_ref <- 21128396
n_het <- 648538
n_alt <- 1001
n_miss <- 51563345
n_all <- n_ref + n_het + n_alt + n_miss # 73,381,805
n_called <- n_ref + n_het + n_alt # 21,818,460
n_nonref <- n_het + n_alt #    649,540
pct <- function(x, tot) 100 * x / tot

## cumulative tops (percent) for each bar's stack
# bar 1 "All": REF, HET, ALT, missing (called block at the bottom)
a_ref <- pct(n_ref, n_all) # 28.85
a_het <- a_ref + pct(n_het, n_all) # 29.73
a_alt <- a_het + pct(n_alt, n_all) # 29.73
called_top <- a_alt # top of the called block
# bar 2 "Called": REF, HET, ALT
c_ref <- pct(n_ref, n_called) # 97.02
c_het <- c_ref + pct(n_het, n_called) # ~100
nonref_bot <- c_ref # bottom of the non-REF block
# bar 3 "Non-REF": HET, ALT
r_het <- pct(n_het, n_nonref) # 99.85

w <- 0.28
rects <- rbind(
  data.frame(x = 1, ymin = 0, ymax = a_ref, fill = GOLD),
  data.frame(x = 1, ymin = a_ref, ymax = a_het, fill = GREEN),
  data.frame(x = 1, ymin = a_het, ymax = a_alt, fill = PURPLE),
  data.frame(x = 1, ymin = called_top, ymax = 100, fill = GREY),
  data.frame(x = 2, ymin = 0, ymax = c_ref, fill = GOLD),
  data.frame(x = 2, ymin = c_ref, ymax = c_het, fill = GREEN),
  data.frame(x = 2, ymin = c_het, ymax = 100, fill = PURPLE),
  data.frame(x = 3, ymin = 0, ymax = r_het, fill = GREEN),
  data.frame(x = 3, ymin = r_het, ymax = 100, fill = PURPLE)
)
rects$xmin <- rects$x - w
rects$xmax <- rects$x + w

## grey "zoom lens" trapezoids: sub-region of the left bar -> full right bar
conn <- rbind(
  data.frame(
    g = "A", x = c(1 + w, 2 - w, 2 - w, 1 + w),
    y = c(0, 0, 100, called_top)
  ),
  data.frame(
    g = "B", x = c(2 + w, 3 - w, 3 - w, 2 + w),
    y = c(nonref_bot, 0, 100, 100)
  )
)

lab <- function(x, y, text, col = INK, size = 4, face = "plain", hj = 0.5) {
  annotate("text",
    x = x, y = y, label = text, color = col, size = size,
    fontface = face, hjust = hj, lineheight = 0.95
  )
}

cc <- function(x) formatC(x, format = "d", big.mark = ",") # 648,538
mk <- function(x) {
  ifelse(x >= 1e6, sprintf("%.1fM", x / 1e6), # 73.3M / 649.5K
    ifelse(x >= 1e3, sprintf("%.1fK", x / 1e3), as.character(x))
  )
}

p <- ggplot() +
  geom_polygon(data = conn, aes(x, y, group = g), fill = GREY, alpha = 0.28) +
  geom_rect(data = rects, aes(
    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
    fill = fill
  ), color = "white", linewidth = 0.4) +
  scale_fill_identity() +
  # bar 1 in-segment labels
  lab(1, a_ref / 2, sprintf("REF-hom\n%.2f%%", a_ref)) +
  lab(1, (called_top + 100) / 2, sprintf("missing\n%.2f%%", 100 - called_top)) +
  # bar 2
  lab(2, c_ref / 2, sprintf("REF-hom\n%.2f%%", c_ref), col = "white", face = "bold") +
  # bar 3
  lab(3, r_het / 2, sprintf("HET\n%.2f%%\n%s", r_het, cc(n_het)), col = "white", face = "bold") +
  # callouts for the invisibly-thin slivers
  annotate("segment",
    x = 3 + w, xend = 3 + w + 0.28, y = (r_het + 100) / 2, yend = 88,
    color = RED, linewidth = 0.5
  ) +
  lab(3 + w + 0.30, 84, sprintf("ALT-hom\n%.2f%%\n%s", 100 - r_het, cc(n_alt)),
    col = RED, size = 3.6, hj = 0
  ) +
  lab(2, c_ref - 5, sprintf("→ non-REF %.2f%%", 100 - nonref_bot), col = INK, size = 3.2) +
  # column captions
  lab(1, -4, "All calls", face = "bold", size = 4.2) +
  lab(1, -9, sprintf("%% of all · %s", mk(n_all)), col = GREY, size = 3.3) +
  lab(2, -4, "Called", face = "bold", size = 4.2) +
  lab(2, -9, sprintf("%% of called · %s", mk(n_called)), col = GREY, size = 3.3) +
  lab(3, -4, "Non-REF called", face = "bold", size = 4.2) +
  lab(3, -9, sprintf("%% of non-REF · %s", mk(n_nonref)), col = GREY, size = 3.3) +
  labs(
    title = "Skim SNP50K genotype calls: almost no donor-hom",
    subtitle = "70% of sites are missing; of the calls 97% are REF and HET is called ~648x more often than ALT-hom."
  ) +
  coord_cartesian(ylim = c(-11, 101), xlim = c(0.55, 4.05), clip = "off") +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, color = INK),
    plot.subtitle = element_text(size = 10.5, color = "grey30"),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(12, 14, 10, 14)
  )

out <- file.path(FIGDIR, "genotype_call_composition.png")
ggsave(out, p, width = 9, height = 5.6, dpi = 200, bg = "white")
cat(sprintf(
  "%% of all:    REF %.2f  HET %.2f  ALT %.4f  miss %.2f\n",
  a_ref, a_het - a_ref, a_alt - a_het, 100 - called_top
))
cat(sprintf(
  "%% of called: REF %.2f  HET %.2f  ALT %.4f\n",
  c_ref, c_het - c_ref, 100 - c_het
))
cat(sprintf(
  "%% non-REF:   HET %.2f  ALT %.2f  (HET:ALT = %.0f:1)\n",
  r_het, 100 - r_het, n_het / n_alt
))
cat("wrote:", out, "\n")
