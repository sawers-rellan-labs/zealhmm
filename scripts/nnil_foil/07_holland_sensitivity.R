#!/usr/bin/env Rscript
# Calibration foil, step 7: WHICH parameter does Holland's model actually care
# about? Uses Holland's OWN full-factorial grid search (File_S04: nir x germ x
# gert x p x r, scored by GBS-vs-chip calls mismatch) to decompose the variance of the
# mismatch objective across the five parameters. Tests the claim that the caller
# is insensitive to r (the recombination/segment-length knob) and driven instead
# by the emission parameters (esp. nir).
#
#   Rscript scripts/nnil_foil/07_holland_sensitivity.R
# Output: data/nnil_foil/holland_param_sensitivity.csv (variance share + marginal
#         range per parameter) + agent/nnil_foil_holland_sensitivity.png

suppressMessages({
  library(data.table)
  library(ggplot2)
})
root <- here::here()
F4 <- file.path(root, "agent/nNIL/File_S04.nNIL_gbs_vs_chip_data_HMMgridSearch.csv")
d <- fread(F4)
params <- c("nir", "germ", "gert", "p", "r")
cat(sprintf(
  "Holland grid: %d parameter combinations (full factorial %s)\n",
  nrow(d), paste(sapply(params, function(p) uniqueN(d[[p]])), collapse = "x")
))

# --- variance decomposition of the mismatch objective across parameters ------
# Full-factorial + balanced -> main-effects ANOVA sums of squares partition the
# objective's variance cleanly across the five knobs. Report each knob's share.
dd <- copy(d)
dd[, (params) := lapply(.SD, factor), .SDcols = params]
fit <- aov(mismatchMean ~ nir + germ + gert + p + r, data = dd)
tab <- as.data.table(summary(fit)[[1]], keep.rownames = "term")
tab[, term := trimws(term)]
ss_tot <- sum(tab$`Sum Sq`)
tab[, var_share := `Sum Sq` / ss_tot]

# --- marginal effect: how much does mean mismatch move across each knob? ------
marg <- rbindlist(lapply(params, function(pp) {
  m <- d[, .(mm = mean(mismatchMean)), by = c(pp)]
  data.table(
    param = pp, n_levels = nrow(m),
    lo = min(d[[pp]]), hi = max(d[[pp]]),
    mismatch_min = min(m$mm), mismatch_max = max(m$mm),
    marginal_range = max(m$mm) - min(m$mm)
  )
}))
out <- merge(
  marg,
  tab[term %in% params, .(param = term, var_share)],
  by = "param"
)[order(-var_share)]
out[, rel_to_best_r := marginal_range / marginal_range[param == "r"]]
fwrite(out, file.path(root, "data/nnil_foil/holland_param_sensitivity.csv"))
cat("\n=== Holland-grid sensitivity of GBS-vs-chip calls mismatch ===\n")
print(out[, .(param,
  levels_lo_hi = sprintf("%.4g..%.4g", lo, hi),
  var_share = round(var_share, 4), marginal_range = signif(marginal_range, 3),
  x_vs_r = round(rel_to_best_r, 1)
)])

# --- figure: marginal mean mismatch vs each parameter ------------------------
plt <- rbindlist(lapply(params, function(pp) {
  m <- d[, .(mm = mean(mismatchMean)), by = c(pp)]
  setnames(m, pp, "level")
  m[, param := pp][]
}))
plt[, param := factor(param,
  levels = out$param,
  labels = sprintf("%s (%.0f%% var)", out$param, 100 * out$var_share)
)]
p_fig <- ggplot(plt, aes(level, mm)) +
  geom_line(linewidth = 0.5, colour = "grey30") +
  geom_point(size = 1.3, colour = "#D55E00") +
  facet_wrap(~param, scales = "free_x", nrow = 1) +
  scale_x_log10() +
  labs(
    x = "parameter value (log scale)", y = "mean GBS-vs-chip calls mismatch",
    title = "Holland's grid: mismatch is driven by nir, not r"
  ) +
  theme_classic(base_size = 9) +
  theme(
    text = element_text(family = "sans"), plot.title = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
OUT <- file.path(root, "agent/nnil_foil_holland_sensitivity.png")
ggsave(OUT, p_fig, width = 190, height = 55, units = "mm", dpi = 300)
cat(sprintf("\nwrote %s\n", OUT))
