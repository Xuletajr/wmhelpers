#' Plot Histogram of MCMC Diagnostics
#'
#' Produces a faceted histogram for a chosen MCMC diagnostic metric
#' (e.g. R-hat, ESS) across model parameter groups (e.g. beta, alpha).
#'
#' @author Thiago Bicudo
#'
#' @param out_bestmodel Fitted model object containing diagnostic lists
#'   (e.g. \code{out_bestmodel$rhat}, \code{out_bestmodel$ESS}).
#'
#' @param metric Character. Name of the diagnostic metric to extract
#'   (e.g. \code{"rhat"}, \code{"ESS"}).
#'
#' @param params Character vector of parameter groups to include
#'   (e.g. \code{c("beta", "alpha")}).
#'
#' @param param_labels Optional named vector for facet labels.
#'   Example: \code{c(beta = "Occupancy", alpha = "Detection")}.
#'
#' @param x_label Character. Label for x-axis.
#'
#' @param threshold Numeric or \code{NULL}. If provided, adds a vertical line.
#'
#' @param output_path Character. If provided, saves the plot as PNG.
#'
#' @param color_fill,color_text,color_grid Character. Plot styling.
#'
#' @param width,height,dpi Numeric. Output dimensions (inches) and resolution.
#'
#' @return Invisibly returns the \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_histogram facet_wrap theme_minimal labs theme element_text element_line ggsave geom_vline
#'
#' @export
plot_mcmc_diagnostic <- function(
    out_bestmodel,
    metric = "rhat",
    params = c("beta", "alpha"),
    param_labels = NULL,
    x_label = NULL,
    threshold = NULL,
    output_path = NULL,
    color_fill = "#CCCCCC",
    color_text = "#1B3A5C",
    color_grid = "#E8EEF4",
    width = 12,
    height = 9,
    dpi = 300
) {

  # ── Validate metric ─────────────────────────────────────────────────────────
  stopifnot(metric %in% names(out_bestmodel))

  metric_list <- out_bestmodel[[metric]]

  stopifnot(all(params %in% names(metric_list)))

  # ── Extract data ────────────────────────────────────────────────────────────
  values <- unlist(metric_list[params])
  group  <- rep(params, times = sapply(metric_list[params], length))

  # ── Labels ──────────────────────────────────────────────────────────────────
  if (!is.null(param_labels)) {
    group <- param_labels[group]
  }

  if (is.null(x_label)) {
    x_label <- metric
  }

  df <- data.frame(
    value = values,
    Parameter = group
  )

  # ── Plot ────────────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(df, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(
      bins = 30,
      fill = color_fill,
      color = color_fill
    ) +
    ggplot2::facet_wrap(~Parameter, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = x_label,
      y = "Frequency"
    ) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = color_text),
      axis.title = ggplot2::element_text(color = color_text),
      panel.grid = ggplot2::element_line(color = color_grid),
      strip.text = ggplot2::element_text(color = color_text)
    )

  # ── Threshold line ──────────────────────────────────────────────────────────
  if (!is.null(threshold)) {
    p <- p + ggplot2::geom_vline(
      xintercept = threshold,
      linetype = "dashed"
    )
  }

  # ── Save ────────────────────────────────────────────────────────────────────
  if (!is.null(output_path)) {
    ggplot2::ggsave(
      filename = output_path,
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
  }

  print(p)

  return(invisible(p))
}
