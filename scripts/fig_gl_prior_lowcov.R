#!/usr/bin/env Rscript
# Figure: "gl_prior_lowcov" — at ~0.4x coverage the genotype PRIOR, not the read,
# decides the single-read call. Four priors compared on one ALT read.
# Spec + calculations: agent/gl_prior_lowcov_doc.md; extends the 2-prior
# scripts/fig_one_alt_read_calls.R to four priors, each tagged with the caller
# that implements it.
#
#   1 ML (flat)            -> GATK single-sample (HaplotypeCaller, argmax PL)
#   2 HWE @ p_hat = 0.125  -> GATK joint (GenotypeGVCFs); p_hat is the THEORETICAL
#                             BC2S3 donor-allele freq f_BB + f_AB/2, computed here
#   3 HWE @ p = 0.015      -> bcftools mpileup | call -m; p is the EMPIRICAL
#                             GL-based cohort mean (a supplied constant)
#   4 BC2S3 design         -> the single-locus NIL genotype freqs (55/64,2/64,7/64);
#                             equals inbreeding-HWE (F=0.857, p=0.125), applicable
#                             via ANGSD -doPost with a per-individual F (ngsF)
#
# Two panels, one column:
#   (a) posterior P(G | one ALT read) as grouped bars, three genotypes per prior
#   (b) P(donor-hom) vs number of ALT reads (r=0), one line per prior
#
# Palette: R color names gold / forestgreen / purple4 (per Fausto's request).

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(here)
})

# Sourced from analysis/genotype_likelihoods_and_hmm.qmd for the `fig` object;
# run standalone (`Rscript scripts/fig_gl_prior_lowcov.R`) to also write the PNG
# and print the tables (guarded on sys.nframe() below).

## ---- palette ---------------------------------------------------------------
GCOL <- c(
  AA = "#FFD700", # gold        REF-hom (B73)
  AB = "#228B22", # forestgreen het
  BB = "#551A8B"
) # purple4     donor-hom (teosinte)
INK <- "#1A1A1A"
GREY <- "#888888"
GENO_LAB <- c(AA = "REF-hom", AB = "het", BB = "donor-hom")

## ---- model -----------------------------------------------------------------
eps <- 0.01

# P(one read | G): is_alt = the observed base is the ALT allele B
read_like <- function(is_alt) {
  if (is_alt) {
    c(AA = eps / 3, AB = 0.5 * (eps / 3) + 0.5 * (1 - eps), BB = 1 - eps)
  } else {
    c(AA = 1 - eps, AB = 0.5 * (1 - eps) + 0.5 * (eps / 3), BB = eps / 3)
  }
}
pileup <- function(a, r) { # a ALT reads, r REF reads
  o <- c(AA = 1, AB = 1, BB = 1)
  if (a) for (i in seq_len(a)) o <- o * read_like(TRUE)
  if (r) for (i in seq_len(r)) o <- o * read_like(FALSE)
  o
}
posterior <- function(a, r, prior) {
  po <- pileup(a, r) * prior
  po / sum(po)
}

## ---- priors ----------------------------------------------------------------
bc2s3 <- c(AA = 55 / 64, AB = 2 / 64, BB = 7 / 64) # 0.859375, 0.03125, 0.109375
p_bc2s3 <- unname(bc2s3["BB"] + bc2s3["AB"] / 2) # THEORETICAL BC2S3 donor freq = 0.125
p_emp <- 0.015 # EMPIRICAL GL-based cohort mean
hwe <- function(p) c(AA = (1 - p)^2, AB = 2 * p * (1 - p), BB = p^2)

priors <- list(
  `1` = c(AA = 1 / 3, AB = 1 / 3, BB = 1 / 3),
  `2` = hwe(p_bc2s3),
  `3` = hwe(p_emp),
  `4` = bc2s3
)

# two-/three-line x labels: prior on top, source, caller (as agreed with Fausto)
xlab <- c(
  `1` = "ML\nGATK single-sample",
  `2` = "HWE @ p̂=0.125\nGATK GenotypeGVCFs",
  `3` = "HWE @ p=0.015\n(GL-based AF)\nbcftools mpileup",
  `4` = "BC2S3 design\n(55/64, 2/64, 7/64)\nANGSD -doPost +F"
)
# short labels for the panel-(b) legend
leglab <- c(
  `1` = "ML flat (GATK 1-sample)",
  `2` = "GATK joint (HWE p̂=0.125)",
  `3` = "bcftools (HWE p=0.015)",
  `4` = "BC2S3 / ANGSD +F"
)
LINECOL <- c(`1` = "#616161", `2` = "#1B7837", `3` = "#4575B4", `4` = "#551A8B")

## ---- panel (a): one-ALT-read posteriors ------------------------------------
dfa <- do.call(rbind, lapply(names(priors), function(k) {
  po <- posterior(1, 0, priors[[k]])
  data.frame(prior = k, g = names(po), post = as.numeric(po))
}))
dfa$prior <- factor(dfa$prior, levels = names(priors))
dfa$g <- factor(dfa$g, levels = c("AA", "AB", "BB"))

# winning-call label above the tallest bar of each group
win <- do.call(rbind, lapply(split(dfa, dfa$prior), function(d) d[which.max(d$post), ]))
win$lab <- GENO_LAB[as.character(win$g)]

dodge <- position_dodge(width = 0.8)
pa <- ggplot(dfa, aes(prior, post, fill = g)) +
  geom_col(position = dodge, width = 0.76, color = "white", linewidth = 0.3) +
  geom_text(
    data = win, aes(label = lab, color = g), position = dodge,
    vjust = -0.5, size = 3.0, fontface = "bold", show.legend = FALSE
  ) +
  scale_fill_manual(values = GCOL, labels = GENO_LAB, name = NULL) +
  scale_color_manual(values = GCOL, guide = "none") +
  scale_x_discrete(labels = xlab) +
  scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  labs(
    title = "(a) One ALT read: the prior decides the call",
    subtitle = "same read, same likelihood PL=[25,3,0]; posterior P(G | one ALT read)",
    x = NULL, y = "posterior  P(G | one ALT read)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 7, lineheight = 0.95, color = INK),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.5, color = "grey30"),
    legend.position = "top", legend.justification = "left",
    legend.key.size = unit(3.4, "mm"), legend.text = element_text(size = 8),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(6, 10, 2, 6)
  )

## ---- panel (b): ALT reads needed to reach donor-hom ------------------------
nn <- 1:12
dfb <- do.call(rbind, lapply(names(priors), function(k) {
  data.frame(
    prior = k, n = nn,
    pbb = vapply(nn, function(a) posterior(a, 0, priors[[k]])["BB"], numeric(1))
  )
}))
dfb$prior <- factor(dfb$prior, levels = names(priors))
# minimum ALT reads to cross P(donor-hom) > 0.5, per prior (for the console table)
mincross <- vapply(names(priors), function(k) {
  ok <- which(dfb$pbb[dfb$prior == k] > 0.5)
  if (length(ok)) nn[min(ok)] else NA_integer_
}, integer(1))

pb <- ggplot(dfb, aes(n, pbb, color = prior)) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = GREY, linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  geom_point(aes(shape = prior), size = 1.9) +
  scale_color_manual(values = LINECOL, labels = leglab, name = "prior (caller)") +
  scale_shape_manual(values = c(16, 17, 15, 18), labels = leglab, name = "prior (caller)") +
  scale_x_continuous(breaks = nn, limits = c(1, 12), expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
  labs(
    title = "(b) ALT reads needed to reach a donor-hom call",
    subtitle = "P(donor-hom) vs ALT read count (r=0); dotted line = 0.5",
    x = "number of ALT reads at the site  (r = 0 REF reads)",
    y = "posterior  P(donor-hom)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.5, color = "grey30"),
    legend.position = "right", legend.key.size = unit(4, "mm"),
    legend.title = element_text(size = 8.5), legend.text = element_text(size = 8),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(6, 10, 4, 6)
  )

fig <- pa / pb + plot_layout(heights = c(1, 0.92))

## ---- write + emit tables (only when run standalone, not when sourced) -------
if (sys.nframe() == 0L) {
  OUTDIR <- here::here("results")
  dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
  OUT <- file.path(OUTDIR, "gl_prior_lowcov.png")
  ggsave(OUT, fig, width = 6.8, height = 8.0, dpi = 200, bg = "white")

  cat("\nOne-ALT-read posteriors (a=1, r=0):\n")
  for (k in names(priors)) {
    po <- posterior(1, 0, priors[[k]])
    g <- names(po)[which.max(po)]
    cat(sprintf(
      "  prior %s  P=(AA %.3f, AB %.3f, BB %.3f)  call=%-9s  AB:BB=%.1f\n",
      k, po["AA"], po["AB"], po["BB"], GENO_LAB[g], po["AB"] / po["BB"]
    ))
  }
  cat("\nMin ALT reads to P(donor-hom) > 0.5:\n")
  for (k in names(priors)) cat(sprintf("  prior %s  %s\n", k, mincross[k]))
  cat(sprintf("\np_hat (BC2S3, computed) = %.4f ; p_emp = %.4f\n", p_bc2s3, p_emp))
  cat("wrote:", OUT, "\n")
}
