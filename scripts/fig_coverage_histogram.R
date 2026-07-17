#!/usr/bin/env Rscript
# Per-sample coverage histogram from missing-data-floor-model.qmd section 3,
# re-rendered with base_size = 25 and exported at 500 x 500 px.
# Data: data/missing_data/wideseq_per_sample.tsv (BZea Wideseq per-sample lambda).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})
theme_cls <- if (requireNamespace("ggpubr", quietly = TRUE)) ggpubr::theme_classic2 else theme_classic

per_sample <- read_tsv("data/missing_data/wideseq_per_sample.tsv", show_col_types = FALSE)
lambda_bar <- sum(per_sample$DEPTH_SUM) / sum(per_sample$VARIANT_COUNT) # ~0.59

p <- per_sample %>%
  ggplot(aes(x = lambda)) +
  xlab("Coverage") +
  ylab("Sample Count") +
  geom_histogram(fill = "grey65", color = "grey30", linewidth = 0.2) +
  geom_vline(xintercept = lambda_bar, color = "red", linewidth = 1) +
  theme_cls(base_size = 25) +
  theme(legend.position = "none")

out <- "nilhmm-paper/figures/coverage_histogram_500.png"
ggsave(out, p, width = 5, height = 5, units = "in", dpi = 100, bg = "white") # 500x500 px
cat(sprintf("lambda_bar = %.3f\nwrote: %s (500x500 px)\n", lambda_bar, out))
