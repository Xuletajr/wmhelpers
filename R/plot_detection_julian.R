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
#' @param start_date Date or character (coercible to \code{Date}). Lower
#'   x-axis limit (the earliest date to display).
#'
#' @param end_date Date or character. Upper x-axis limit (the latest date to
#'   display).
#'
#' @param origin_date Date or character. The calendar date that \code{seq_julian
#'   == 1} corresponds to -- i.e. the arithmetic anchor used to convert raw
#'   Julian-day values into calendar dates. This is \strong{not} always the
#'   same as \code{start_date}; see Details. Defaults to \code{start_date} for
#'   backward compatibility with the year-boundary use case described below.
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
#' @details
#' \code{pred_det_df$seq_julian} (produced by \code{predict_detection_julian()})
#' holds raw Julian-day values, and the date shown for each point is computed
#' as \code{origin_date + (seq_julian - 1)}. What \code{origin_date} should be
#' depends on how \code{julian_day_raw} was built upstream, in the project's
#' own data-preparation step -- there are two legitimate conventions:
#' \itemize{
#'   \item \strong{Calendar day-of-year, no year wraparound} (e.g. within-year
#'     seasons that don't cross a year boundary): \code{seq_julian} is the
#'     ordinary day-of-year (1-366), and \code{origin_date} must be
#'     \strong{January 1st} of a reference year -- \emph{not} \code{start_date}.
#'     Using \code{start_date} here double-counts the offset, since
#'     \code{start_date} is itself already Jan 1 plus some number of days.
#'   \item \strong{Project spans a year boundary} (e.g. a Nov-Mar season) and
#'     \code{julian_day_raw} was deliberately re-anchored upstream so that day
#'     1 = the project's start date (avoiding the 365 -> 0 -> 90 jump a
#'     calendar day-of-year would otherwise produce): \code{origin_date} is
#'     the actual project start date, i.e. the same value as \code{start_date}.
#' }
#' \code{start_date} and \code{end_date} always control only the x-axis
#' display limits (via \code{scale_x_date(limits = ...)}) and are independent
#' of \code{origin_date}; conflating the two (as this function used to, before
#' \code{origin_date} was a separate argument) silently double-counts the
#' offset whenever the first convention above applies.
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line theme_minimal labs ggsave scale_x_date
#'
#' @export
plot_detection_julian <- function(
    pred_det_df,
    start_date,
    end_date,
    origin_date = start_date,
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

  start_date  <- as.Date(start_date)
  end_date    <- as.Date(end_date)
  origin_date <- as.Date(origin_date)

  if (is.na(start_date) || is.na(end_date) || is.na(origin_date)) {
    stop("start_date, end_date, and origin_date must be valid Date or coercible to Date")
  }

  # ── Convert Julian → Date ───────────────────────────────────────────────────
  # origin_date (NOT start_date) is the arithmetic anchor -- see Details for
  # which value to pass depending on how seq_julian was constructed upstream.
  pred_det_df$date <- origin_date + (pred_det_df$seq_julian - 1)

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
