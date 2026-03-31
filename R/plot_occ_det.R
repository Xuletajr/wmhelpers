#' Plot Mean Occupancy vs Mean Detection Probability
#'
#' Scatter plot of species-level mean occupancy (x-axis) against mean detection
#' probability (y-axis), including a linear regression line and annotation with
#' R² and p-value.
#'
#' @author Thiago Bicudo
#'
#' @param export_occ_prob A tidy data frame from \code{extract_prob} containing
#'   occupancy estimates. Must include columns \code{species} and \code{mean}.
#'
#' @param export_det_prob A tidy data frame from \code{extract_prob} containing
#'   detection estimates. Must include columns \code{species} and \code{mean}.
#'
#' @param color_line Character. Color of the regression line.
#'   Default is \code{"#1B3A5C"}.
#'
#' @param output_path Character. If provided, saves the plot as a PNG file.
#'   Default is \code{NULL}.
#'
#' @param width,height,dpi Numeric. Dimensions (in inches) and resolution for
#'   saved PNG. Defaults: \code{width = 9}, \code{height = 9}, \code{dpi = 600}.
#'
#' @return Invisibly returns the fitted \code{lm} object.
#'
#' @importFrom dplyr inner_join rename
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth theme_minimal labs annotate
#'
#'#' @author Thiago Bicudo
#'
#' @export
plot_occ_vs_det <- function(
    export_occ_prob,
    export_det_prob,
    color_line = "#1B3A5C",
    output_path = NULL,
    width = 9,
    height = 9,
    dpi = 600
) {

  # ── Input validation ────────────────────────────────────────────────────────
  stopifnot(
    is.data.frame(export_occ_prob),
    is.data.frame(export_det_prob),
    all(c("species", "mean") %in% names(export_occ_prob)),
    all(c("species", "mean") %in% names(export_det_prob))
  )

  # ── Prepare data ────────────────────────────────────────────────────────────
  df <- dplyr::inner_join(
    export_occ_prob,
    export_det_prob,
    by = "species"
  ) %>%
    dplyr::rename(
      occ = mean.x,
      det = mean.y
    )

  # ── Fit linear model ────────────────────────────────────────────────────────
  mod <- lm(det ~ occ, data = df)

  r2 <- summary(mod)$r.squared
  pval <- summary(mod)$coefficients[2, 4]

  label_text <- paste0(
    "R² = ", round(r2, 3),
    "\np = ", signif(pval, 3)
  )

  # ── Plot ────────────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(df, ggplot2::aes(x = occ, y = det)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = color_line) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Mean Occupancy Probability",
      y = "Mean Detection Probability"
    ) +
    ggplot2::annotate(
      "text",
      x = min(df$occ, na.rm = TRUE),
      y = max(df$det, na.rm = TRUE),
      label = label_text,
      hjust = 0,
      vjust = 1
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

  return(invisible(mod))
}
