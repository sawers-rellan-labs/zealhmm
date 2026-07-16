# Per-sample missing-data model: the exponential-with-floor fit and its two
# companion plots (coverage histogram + observed-vs-expected validation panel),
# factored out of the missing-data notebooks so the fit and figure that produced
# the ZEAL slide deck (data/20250819_Fausto.pdf, p.2) live in one place and never
# drift from a hand-retyped copy.
#
# Ported verbatim from the original modeling_functions.R. Assumes the caller has
# already loaded: dplyr/magrittr (%>%), ggplot2, minpack.lm (nlsLM), ggpubr,
# cowplot. viridis scales come from ggplot2.

#' Fit the exponential-with-floor missingness model to a per-sample data frame
#'
#' Expects columns `lambda` (mean coverage) and `missing_obs` (observed
#' missingness). Fits `missing_obs ~ pi + (1 - pi) * exp(-k * lambda)`.
#'
#' @param per_sample Data frame with `lambda` and `missing_obs`.
#' @return list(fit, params = list(pi, k), data = per_sample + `missing_exp_floor`).
#' @export
fit_exp_floor_model <- function(per_sample) {
  per_sample <- per_sample %>%
    mutate(missing_poisson = exp(-lambda))

  dat <- per_sample %>%
    filter(
      is.finite(lambda), is.finite(missing_obs),
      missing_obs > 0, missing_obs < 1
    )

  # Starting values: linear fits for the floor intercept and the decay rate
  pi_est <- coef(lm(missing_obs ~ missing_poisson, data = per_sample))[1]
  k_est <- -coef(lm(log(missing_obs) ~ lambda, data = dat))[2]

  fit <- nlsLM(
    missing_obs ~ pi + (1 - pi) * exp(-k * lambda),
    data  = dat,
    start = list(pi = unname(pi_est), k = unname(k_est)),
    lower = c(0, 0),
    upper = c(1, 10)
  )

  params <- as.list(coef(fit))
  per_sample$missing_exp_floor <- predict(fit, newdata = per_sample)

  list(fit = fit, params = params, data = per_sample)
}

#' Coverage histogram for a per-sample data frame (used as an inlay)
#'
#' @param per_sample Data frame with a `lambda` (coverage) column.
#' @return A ggplot object; a dashed red line marks the mean coverage.
#' @export
generate_coverage_histogram <- function(per_sample) {
  mean_lambda <- mean(per_sample$lambda, na.rm = TRUE)

  per_sample %>%
    ggplot(aes(x = lambda)) +
    geom_histogram(aes(fill = after_stat(x))) +
    geom_vline(xintercept = mean_lambda, color = "red", linetype = "dashed") +
    scale_fill_viridis_c() +
    scale_x_continuous(limits = c(0, 1.5)) +
    labs(x = "Coverage", y = "Count") +
    ggpubr::theme_classic2(base_size = 8) +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7)
    )
}

#' Observed-vs-expected missingness validation panel
#'
#' Scatter of observed vs floor-model-expected missingness coloured by coverage,
#' with the 1:1 line (red), the fitted floor `pi` (dashed grey vertical), the
#' mean observed missingness (dashed red horizontal), and the coverage histogram
#' as a lower-right inlay. Title carries the fitted pi/k and mean missingness.
#'
#' @param model_results Output of `fit_exp_floor_model()`.
#' @param coverage_hist_plot A ggplot from `generate_coverage_histogram()`.
#' @param plot_title Panel title (dataset name).
#' @param xy_limits Shared x/y axis limits.
#' @param legend_breaks,legend_labels Coverage colourbar breaks/labels.
#' @return A cowplot grob combining the scatter and the histogram inlay.
#' @export
generate_floor_model_plot <- function(model_results, coverage_hist_plot, plot_title,
                                      xy_limits = c(0.1, 1),
                                      legend_breaks = c(0, 0.5, 1, 1.5, 2),
                                      legend_labels = as.character(legend_breaks)) {
  p_scatter_floor <- model_results$data %>%
    ggplot(aes(x = missing_exp_floor, y = missing_obs, col = lambda)) +
    geom_point(alpha = 0.8) +
    geom_vline(
      xintercept = model_results$params$pi,
      linetype = "dashed", color = "gray50"
    ) +
    geom_hline(
      yintercept = mean(model_results$data$missing_obs),
      linetype = "dashed", color = "red"
    ) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    scale_color_viridis_c(
      limits = range(legend_breaks),
      breaks = legend_breaks,
      labels = legend_labels
    ) +
    guides(col = guide_colorbar(title = "Coverage")) +
    coord_fixed(ratio = 1) +
    scale_x_continuous(limits = xy_limits) +
    scale_y_continuous(limits = xy_limits) +
    labs(
      title = plot_title,
      subtitle = paste0(
        "π = ", round(model_results$params$pi, 3),
        ", k = ", round(model_results$params$k, 3),
        ", missing =", round(mean(model_results$data$missing_obs, na.rm = TRUE), 3)
      ),
      x = "Expected Missing (Floor Model)",
      y = "Observed Missing"
    ) +
    ggpubr::theme_classic2(base_size = 20) +
    theme(legend.position = "right")

  cowplot::ggdraw() +
    cowplot::draw_plot(p_scatter_floor) +
    cowplot::draw_plot(cowplot::as_grob(coverage_hist_plot),
      x = 0.45, y = 0.23, width = 0.25, height = 0.25
    )
}
