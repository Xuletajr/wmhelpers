#' Plot Histogram of R-hat Values
#'
#' Produces a faceted histogram of R-hat values for beta (occupancy) and
#' alpha (detection) parameters and optionally exports it to PNG.
#'
#' @author Thiago Bicudo
#'
#' @param out_bestmodel Fitted model object with \code{rhat$beta} and
#'   \code{rhat$alpha} numeric vectors.
#'
#' @param output_path Character. Full path including filename for PNG export.
#'   If \code{NULL} (default), the plot is printed to the active device only.
#'
#' @param color_fill Character. Fill and border colour for histogram bars.
#'   Defaults to \code{"#CCCCCC"}.
#'
#' @param color_text Character. Axis text colour. Defaults to \code{"#1B3A5C"}.
#'
#' @param color_grid Character. Panel grid colour. Defaults to \code{"#E8EEF4"}.
#'
#' @param width,height,dpi Numeric. PNG dimensions in inches and resolution.
#'
#' @return Invisibly returns the \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_histogram facet_wrap theme_minimal labs theme element_text element_line ggsave geom_vline
#'
#' @export
plot_rhat_hist <- function(
    out_bestmodel,
    output_path = NULL,
    color_fill = "#CCCCCC",
    color_text = "#1B3A5C",
    color_grid = "#E8EEF4",
    width = 12,
    height = 9,
    dpi = 300
) {

  # ── Input validation ────────────────────────────────────────────────────────
  stopifnot(
    !is.null(out_bestmodel$rhat),
    "beta" %in% names(out_bestmodel$rhat),
    "alpha" %in% names(out_bestmodel$rhat)
  )

  beta_rhat  <- out_bestmodel$rhat$beta
  alpha_rhat <- out_bestmodel$rhat$alpha

  # ── Build tidy data ─────────────────────────────────────────────────────────
  df <- data.frame(
    Rhat = c(beta_rhat, alpha_rhat),
    Parameter = rep(
      c("Occupancy (beta)", "Detection (alpha)"),
      times = c(length(beta_rhat), length(alpha_rhat))
    )
  )

  # ── Plot ────────────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Rhat)) +
    ggplot2::geom_histogram(
      bins = 30,
      fill = color_fill,
      color = color_fill
    ) +
    ggplot2::geom_vline(
      xintercept = 1.1,
      linetype = "dashed"
    ) +
    ggplot2::facet_wrap(~Parameter, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "R-hat",
      y = "Frequency"
    ) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = color_text),
      axis.title = ggplot2::element_text(color = color_text),
      panel.grid = ggplot2::element_line(color = color_grid),
      strip.text = ggplot2::element_text(color = color_text)
    )

  # ── Save if requested ───────────────────────────────────────────────────────
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
