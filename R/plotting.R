# Paper-figure styling: the REF/HET/ALT palette and the chromosome-painting
# track plot, factored out of the zealtiger sanity paint so every note that
# paints ancestry (main + supplementary) shares one look.

#' The canonical REF/HET/ALT fill palette
#' @export
state_palette <- function() {
  c(REF = "gold", HET = "springgreen4", ALT = "purple4")
}

#' Chromosome-painting track plot from common-schema segments
#'
#' Rows = tracks (source x caller, or samples); columns = chromosomes; each
#' segment a colored rectangle. Expects a `rowlab` factor (row order) and a
#' `method`/track label already attached; states integer-coded 0/1/2.
#'
#' @param df Common-schema segments with extra `rowlab` (factor) and `track`
#'   columns. Rows are drawn top-to-bottom in `levels(df$rowlab)` order.
#' @param title,subtitle Plot labels.
#' @param base_size Base font size.
#' @return A ggplot object.
#' @export
paint_tracks <- function(df, title = NULL, subtitle = NULL, base_size = 8) {
  df <- data.table::copy(data.table::as.data.table(df))
  df[, chr := factor(chr, levels = 1:10, labels = paste0("chr", 1:10))]
  df[, state := factor(state, 0:2, c("REF", "HET", "ALT"))]
  ggplot2::ggplot(df, ggplot2::aes(
    xmin = start_bp / 1e6, xmax = end_bp / 1e6,
    ymin = 0.02, ymax = 0.98, fill = state
  )) +
    ggplot2::geom_rect() +
    ggplot2::facet_grid(rowlab ~ chr,
      scales = "free_x", space = "free_x", switch = "y"
    ) +
    ggplot2::scale_fill_manual(values = state_palette(), drop = FALSE) +
    ggplot2::labs(
      title = title, subtitle = subtitle,
      x = "Position (Mb)", fill = "Genotype"
    ) +
    ggpubr::theme_classic2(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.spacing.x = ggplot2::unit(0.04, "lines"),
      panel.spacing.y = ggplot2::unit(0.15, "lines"),
      strip.background = ggplot2::element_blank(),
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, lineheight = 0.9),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
}
