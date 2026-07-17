#!/usr/bin/env Rscript
# Figure: "1 ALT read — genotype calls".
# Reproducible rebuild of nilhmm-paper/figures/one_ALT_read_Calls.png, condensing
# analysis/genotype_likelihoods_and_hmm.qmd sections 1-2 (likelihood -> prior ->
# posterior) for a single ALT read at epsilon = 0.01.
#
# Key change vs the semi-manual original: the HWE prior is computed at the BC2S3
# DESIGN donor-allele frequency p = 0.125 (not the cohort-empirical 0.0152), so
# HWE and the breeding genotype frequencies are compared at the SAME allele
# frequency. The flip is then purely genotype structure (HWE vs selfing), not
# allele frequency: HWE forces HET ~7:1, breeding forces ALT ~7:1.
# NOTE: this makes the single-read HET:ALT ratio ~7:1, not the ~66:1 that appears
# in the text at p=0.0152 -- update the prose if adopting p=0.125 there too.
#
# Emits two versions:
#   one_ALT_read_Calls_2priors.png  — likelihood x {HWE, breeding} priors -> posteriors
#   one_ALT_read_Calls_HWE.png      — likelihood x HWE prior -> posterior (condensed)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggtext)
})

here <- tryCatch(dirname(sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)
)), error = function(e) ".")
if (length(here) == 0 || here == "") here <- "scripts"
FIGDIR <- normalizePath(file.path(here, "..", "nilhmm-paper", "figures"), mustWork = FALSE)

## ---- palette (matches the genotype paint) --------------------------------
GOLD <- "#E0A81C"
GREEN <- "#2E9B57"
PURPLE <- "#6B3FA0"
INK <- "#1A1A1A"
RED <- "#C0392B"
STATE <- c(AA = GOLD, AB = GREEN, BB = PURPLE)

## ---- model -----------------------------------------------------------------
eps <- 0.01
p <- 0.125 # BC2S3 design donor-allele frequency

# single-ALT-read genotype likelihoods  P(read = B | G)
L <- c(
  AA = eps / 3,
  AB = 0.5 * eps / 3 + 0.5 * (1 - eps),
  BB = 1 - eps
)
L_norm <- L / sum(L) # normalized likelihood (flat prior) for display

# priors
prior_hwe <- c(AA = (1 - p)^2, AB = 2 * p * (1 - p), BB = p^2) # HWE at p=0.125
prior_breed <- c(AA = 0.859375, AB = 0.03125, BB = 0.109375) # BC2S3 gt freqs

posterior <- function(prior) {
  post <- L * prior
  post / sum(post)
}
post_hwe <- posterior(prior_hwe)
post_breed <- posterior(prior_breed)

## ---- one bar panel ---------------------------------------------------------
fmt <- function(v) formatC(v, format = "g", digits = 3)

bar_panel <- function(vals, concept, caller = NULL, highlight = NULL, ymax = 1) {
  df <- data.frame(
    g = factor(names(vals), levels = c("AA", "AB", "BB")),
    v = as.numeric(vals)
  )
  df$lc <- INK
  if (!is.null(highlight)) df$lc[highlight] <- RED # winning genotype label -> red
  # title: concept in ink, caller in red (inline, one line)
  ttl <- if (is.null(caller)) {
    concept
  } else {
    sprintf("%s &middot; <span style='color:%s'>%s</span>", concept, RED, caller)
  }
  ylab_y <- -0.05 # genotype label sits just below the baseline (drawn in-panel)
  ybot <- -0.19 # red box bottom, below the label
  gg <- ggplot(df, aes(g, v, fill = g)) +
    geom_col(width = 0.66, color = "grey25", linewidth = 0.25)
  if (!is.null(highlight)) { # red box encloses frequency + bar + genotype label
    i <- highlight
    gg <- gg + annotate("rect",
      xmin = i - 0.42, xmax = i + 0.42,
      ymin = ybot, ymax = df$v[i] + 0.11,
      fill = NA, color = RED, linewidth = 1.1
    )
  }
  gg +
    geom_text(aes(label = fmt(v)), vjust = -0.5, size = 3.2, color = INK) +
    geom_text(aes(y = ylab_y, label = g, color = lc),
      vjust = 1, size = 3.6,
      fontface = "bold"
    ) +
    scale_color_identity(guide = "none") +
    scale_fill_manual(values = STATE, guide = "none") +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(ybot - 0.02, ymax), clip = "off") +
    labs(title = ttl, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(), axis.ticks = element_blank(),
      # 80% of the "One ALT read" title (7.5 mm ≈ 21 pt): ~17 pt
      plot.title = element_markdown(
        color = INK, face = "bold", size = 17,
        hjust = 0.5, margin = margin(b = 1)
      ),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(2, 6, 2, 6)
    )
}

## ---- panels ----------------------------------------------------------------
p_lik <- bar_panel(L_norm, "Likelihood", caller = "1 sample GATK", highlight = 3)
p_hwe_p <- bar_panel(prior_hwe, "×  HWE prior  (p = 0.125)")
p_brd_p <- bar_panel(prior_breed, "×  Breeding prior  (BC2S3)")
p_hwe_o <- bar_panel(post_hwe, "HWE posterior", caller = "pileup", highlight = 2)
p_brd_o <- bar_panel(post_breed, "Breeding posterior", highlight = 3)

## title text for the empty top-left quarter (version B), no em dash
title_panel <- ggplot() +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void() +
  annotate("text",
    x = 0.03, y = 0.94, label = "One ALT read", hjust = 0, vjust = 1,
    fontface = "bold", size = 7.5, color = INK
  ) +
  annotate("text",
    x = 0.03, y = 0.75, label = "genotype calls", hjust = 0, vjust = 1,
    fontface = "bold", size = 7.5, color = INK
  ) +
  annotate("text",
    x = 0.03, y = 0.44,
    label = "The HWE prior, not the read,\ndecides the call.",
    hjust = 0, vjust = 1, size = 4.3, color = "grey30"
  ) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 8, 8, 12)
  )

annot <- function(pw, sub) {
  pw + plot_annotation(
    title = "One ALT read: genotype calls",
    subtitle = sub,
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, color = INK),
      plot.subtitle = element_text(size = 10, color = "grey30"),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )
}

## version A: both priors (likelihood spans col 1; priors top, posteriors bottom)
verA <- annot(
  (p_lik | (p_hwe_p / p_hwe_o) | (p_brd_p / p_brd_o)) + plot_layout(widths = c(1, 1, 1)),
  "Same donor-allele frequency (0.125): HWE forces HET, the breeding design forces ALT; the flip is genotype structure, not allele frequency."
)

## version B: HWE only, original composite layout with the breeding column dropped.
## 2x2 grid (all three bar panels same size -> identical scale); title fills the
## empty top-left quarter, as in the hand-drawn schema:
##   title      | HWE prior (top-right)
##   likelihood | HWE posterior (bottom-right)
verB <- (title_panel + p_lik + p_hwe_p + p_hwe_o + plot_layout(design = "AC\nBD")) +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", color = NA)))

## version C: the same, for one REF read. Here the prior REINFORCES the read
## rather than flipping it: likelihood -> AA (2:1), HWE posterior -> AA (~7:1).
L_ref <- c(
  AA = 1 - eps, # P(read = A | G)
  AB = 0.5 * (1 - eps) + 0.5 * eps / 3,
  BB = eps / 3
)
L_ref_norm <- L_ref / sum(L_ref)
post_hwe_ref <- {
  pr <- L_ref * prior_hwe
  pr / sum(pr)
}

p_lik_ref <- bar_panel(L_ref_norm, "Likelihood", caller = "1 sample GATK", highlight = 1)
p_hwe_o_ref <- bar_panel(post_hwe_ref, "HWE posterior", caller = "pileup", highlight = 1)

title_panel_ref <- ggplot() +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void() +
  annotate("text",
    x = 0.03, y = 0.94, label = "One REF read", hjust = 0, vjust = 1,
    fontface = "bold", size = 7.5, color = INK
  ) +
  annotate("text",
    x = 0.03, y = 0.75, label = "genotype calls", hjust = 0, vjust = 1,
    fontface = "bold", size = 7.5, color = INK
  ) +
  annotate("text",
    x = 0.03, y = 0.44,
    label = "Likelihood and prior agree:\nREF-hom.",
    hjust = 0, vjust = 1, size = 4.3, color = "grey30"
  ) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 8, 8, 12)
  )

verC <- (title_panel_ref + p_lik_ref + p_hwe_p + p_hwe_o_ref + plot_layout(design = "AC\nBD")) +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", color = NA)))

## ---- write -----------------------------------------------------------------
outA <- file.path(FIGDIR, "one_ALT_read_Calls_2priors.png")
outB <- file.path(FIGDIR, "one_ALT_read_Calls_HWE.png")
outC <- file.path(FIGDIR, "one_REF_read_Calls_HWE.png")
ggsave(outA, verA, width = 11, height = 4.6, dpi = 200, bg = "white")
ggsave(outB, verB, width = 7, height = 7, dpi = 200, bg = "white") # square
ggsave(outC, verC, width = 7, height = 7, dpi = 200, bg = "white") # square

cat(sprintf(
  "ALT read HWE posterior: AA=%.3f AB=%.3f BB=%.3f  (AB:BB = %.1f:1)\n",
  post_hwe["AA"], post_hwe["AB"], post_hwe["BB"], post_hwe["AB"] / post_hwe["BB"]
))
cat(sprintf(
  "REF read HWE posterior: AA=%.3f AB=%.3f BB=%.4f  (AA:AB = %.1f:1)\n",
  post_hwe_ref["AA"], post_hwe_ref["AB"], post_hwe_ref["BB"],
  post_hwe_ref["AA"] / post_hwe_ref["AB"]
))
cat("wrote:\n  ", outA, "\n  ", outB, "\n  ", outC, "\n")
