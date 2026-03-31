#' Plot Community-Level Detection Probability Over Survey Dates
#'
#' Converts a sequence of Julian days to calendar dates and plots the
#' community-level detection probability as a line with a credible interval
#' ribbon.
#'
#' @author Thiago Bicudo
#'
#' @param pred_det_df Data frame as returned by
#'   \code{predict_detection_julian()$pred_det_df}. Must contain columns:
#'   \code{seq_julian}, \code{mean}, \code{criL}, \code{criU}.
#'
#' @param start_date Date or character (coercible to \code{Date}). Reference
#'   start date corresponding to Julian day = 1.
#'
#' @param end_date Date or character. End date for x-axis limits.
#'
#' @param color_ribbon Character. Fill color for the credible interval ribbon.
#'   Default is \code{"#CCCCCC"}.
#'
#' @param color_line Character. Line color. Default is \code{"#1B3A5C"}.
#'
#' @param ribbon_alpha Numeric. Transparency of ribbon. Default is \code{0.55}.
#'
#' @param output_path Character. If provided, saves the plot as PNG.
#'   Default is \code{NULL}.
#'
#' @param width,height Numeric. Plot dimensions in cm.
#'
#' @param dpi Integer. Resolution. Default is \code{300}.
#'
#' @return Invisibly returns the \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line theme_minimal labs ggsave scale_x_date
#'
#' @export
plot_detection_julian <- function(
    pred_det_df,
    start_date,
    end_date,
    color_ribbon = "#CCCCCC",
    color_line = "#1B3A5C",
    ribbon_alpha = 0.55,
    output_path = NULL,
    width = 28,
    height = 24,
    dpi = 300
) {

  # ── Input validation ────────────────────────────────────────────────────────
  stopifnot(
    is.data.frame(pred_det_df),
    all(c("seq_julian", "mean", "criL", "criU") %in% names(pred_det_df))
  )

  start_date <- as.Date(start_date)
  end_date   <- as.Date(end_date)

  if (is.na(start_date) || is.na(end_date)) {
    stop("start_date and end_date must be valid Date or coercible to Date")
  }

  # ── Convert Julian → Date ───────────────────────────────────────────────────
  pred_det_df$date <- start_date + (pred_det_df$seq_julian - 1)

  # ── Plot ────────────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(pred_det_df, ggplot2::aes(x = date, y = mean)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = criL, ymax = criU),
      fill = color_ribbon,
      alpha = ribbon_alpha
    ) +
    ggplot2::geom_line(color = color_line, linewidth = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Date",
      y = "Detection probability"
    ) +
    ggplot2::scale_x_date(limits = c(start_date, end_date))

  # ── Save if requested ───────────────────────────────────────────────────────
  if (!is.null(output_path)) {
    ggplot2::ggsave(
      filename = output_path,
      plot = p,
      width = width / 2.54,   # cm → inches
      height = height / 2.54,
      dpi = dpi
    )
  }

  print(p)

  return(invisible(p))
}
