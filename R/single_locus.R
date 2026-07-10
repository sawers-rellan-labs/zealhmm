# Single-locus validation helpers — ported verbatim from the zealtiger
# exploratory repo (R/07_single_locus.R + truth_segments_from_dosage from
# R/05_make_rtiger_input.R). Used by analysis/single-locus-validation.qmd.
#
# The BC2S3 simulation's single-locus MARGINAL genotype distribution equals the
# breeding-scheme transition-matrix expectation, exactly and independent of the
# map and of interference m (a marginal quantity). These functions compute that
# expectation, tabulate per-NIL Mb-based REF/Het/ALT fractions, and test the
# population means against it with the design-correct Hotelling chi-square
# (sample = unit; EMPIRICAL between-sample covariance, NOT per-marker binomial).

suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

#' Single-locus genotype expectation for a BC(n_bc) S(n_self) scheme
#' @return Named numeric `c(REF, Het, ALT)` summing to 1.
single_locus_p0 <- function(n_bc = 2L, n_self = 3L) {
  B <- rbind(c(0, 1, 0), c(0, 0.5, 0.5), c(0, 0, 1)) # backcross to aa
  S <- rbind(c(1, 0, 0), c(0.25, 0.5, 0.25), c(0, 0, 1)) # selfing
  v <- c(0, 1, 0) # F1 = Aa; order (AA,Aa,aa)
  for (i in seq_len(n_bc)) v <- as.numeric(v %*% B)
  for (i in seq_len(n_self)) v <- as.numeric(v %*% S)
  c(REF = v[3], Het = v[2], ALT = v[1])
}

#' Run-length truth segments from a per-marker dosage vector (0/1/2).
truth_segments_from_dosage <- function(markers, dosage) {
  df <- markers %>%
    dplyr::transmute(.data$chr, .data$bp, state = dosage) %>%
    dplyr::arrange(.data$chr, .data$bp)
  df %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::mutate(run = cumsum(c(TRUE, diff(.data$state) != 0))) %>%
    dplyr::group_by(.data$chr, .data$run, .data$state) %>%
    dplyr::summarise(
      start_bp = min(.data$bp), end_bp = max(.data$bp),
      n_markers = dplyr::n(), .groups = "drop"
    ) %>%
    dplyr::arrange(.data$chr, .data$start_bp) %>%
    dplyr::select("chr", "start_bp", "end_bp", "state", "n_markers")
}

#' Per-NIL Mb-based REF/Het/ALT fractions from a dosage vector.
#' @return Named numeric `c(REF, Het, ALT)` Mb fractions summing to 1.
state_fractions_mb <- function(markers, dosage) {
  seg <- truth_segments_from_dosage(markers, dosage)
  seg$mb <- (seg$end_bp - seg$start_bp) / 1e6
  tot <- tapply(seg$mb, factor(seg$state, levels = 0:2), sum)
  tot[is.na(tot)] <- 0
  f <- as.numeric(tot) / sum(tot)
  c(REF = f[1], Het = f[2], ALT = f[3])
}

#' Design-correct goodness-of-fit of population fractions to an expectation.
#' Sample = unit; Hotelling/Wald chi-square on (Het, ALT) with EMPIRICAL
#' between-sample covariance (df = 2); per-state Wald 95% CI, SE = sd(f_i)/sqrt(N).
#' @return list: T2, df, p_chisq, F_stat, p_F, table (per-state CIs), N.
hotelling_fractions <- function(F, p0, comp = c("Het", "ALT")) {
  stopifnot(ncol(F) == 3, all(c("REF", "Het", "ALT") %in% colnames(F)))
  N <- nrow(F)
  phat <- colMeans(F)
  d <- phat[comp] - p0[comp]
  S2 <- stats::cov(F[, comp])
  T2 <- as.numeric(N * t(d) %*% solve(S2) %*% d)
  Ff <- (N - 2) / (2 * (N - 1)) * T2
  cv <- stats::cov(F)
  se <- sqrt(diag(cv) / N)
  z <- stats::qnorm(0.975)
  st <- c("REF", "Het", "ALT")
  tab <- data.frame(
    state = st, expected = as.numeric(p0[st]), estimate = as.numeric(phat[st]),
    se = as.numeric(se[st]),
    ci_lo = as.numeric(phat[st] - z * se[st]),
    ci_hi = as.numeric(phat[st] + z * se[st]),
    sd_fi = as.numeric(sqrt(diag(cv))[st])
  )
  tab$in_ci <- tab$expected >= tab$ci_lo & tab$expected <= tab$ci_hi
  list(
    T2 = T2, df = 2L, p_chisq = stats::pchisq(T2, 2, lower.tail = FALSE),
    F_stat = Ff, p_F = stats::pf(Ff, 2, N - 2, lower.tail = FALSE),
    table = tab, N = N
  )
}

#' Two-panel forest plot of the fraction validation (A: log magnitude; B: z).
#' @return A patchwork object (A | B).
plot_fraction_forest <- function(tab, title = NULL) {
  ord <- c("REF", "ALT", "Het") # high -> low, top -> bottom
  tab$state <- factor(tab$state, levels = rev(ord))
  tab$z <- (tab$estimate - tab$expected) / tab$se
  tab$pass <- ifelse(tab$in_ci, "in CI", "out")
  pal <- c("in CI" = "#2c7fb8", "out" = "#d7301f")

  pA <- ggplot2::ggplot(tab, ggplot2::aes(y = .data$state)) +
    ggplot2::geom_segment(ggplot2::aes(
      x = .data$ci_lo, xend = .data$ci_hi,
      yend = .data$state, colour = .data$pass
    ), linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(x = .data$estimate, colour = .data$pass), size = 2.6) +
    ggplot2::geom_point(ggplot2::aes(x = .data$expected),
      shape = 124, size = 6,
      colour = "black"
    ) + # expected = vertical tick
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::labs(
      x = "Mb fraction (log scale) — point=estimate ±95% CI, tick=expected",
      y = NULL, subtitle = "A. Magnitude"
    ) +
    ggplot2::theme_bw(base_size = 11)

  pB <- ggplot2::ggplot(tab, ggplot2::aes(y = .data$state)) +
    ggplot2::annotate("rect",
      xmin = -1.96, xmax = 1.96, ymin = -Inf, ymax = Inf,
      fill = "grey85", alpha = 0.7
    ) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::geom_segment(ggplot2::aes(
      x = 0, xend = .data$z, yend = .data$state,
      colour = .data$pass
    ), linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(x = .data$z, colour = .data$pass), size = 2.6) +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::labs(
      x = "standardized deviation (estimate − expected)/SE",
      y = NULL, subtitle = "B. Inside 95% CI?  (band = ±1.96)"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank())

  p <- patchwork::wrap_plots(pA, pB, widths = c(1, 1))
  if (!is.null(title)) p <- p + patchwork::plot_annotation(title = title)
  p
}
