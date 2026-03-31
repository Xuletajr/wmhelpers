#' Build Occupancy Prediction Plot List
#'
#' Produces a named list of \code{ggplot2} objects showing occupancy
#' probability responses across predictor gradients. Each plot corresponds
#' to a combination of variable × year × (optionally) class.
#'
#' @param pred_occ_df Tidy data frame returned by
#'   \code{\link{summarise_occ_predictions}}. Must include columns:
#'   \code{variable}, \code{value}, \code{mean}, \code{criL},
#'   \code{criU}, and optionally \code{year} and \code{class}.
#' @param type_levels Optional character vector of class/type levels
#'   (e.g. habitat types). Defaults to \code{NULL}.
#' @param color_ribbon Character. Ribbon fill colour.
#' @param color_line Character. Line colour.
#' @param color_text Character. Text/axis colour.
#' @param color_strip_bg Character. Facet strip background colour.
#' @param ribbon_alpha Numeric. Transparency of ribbon.
#'   Defaults to \code{0.55}.
#' @param strip_class_prefix Character string to remove from class labels.
#'
#' @return Named list of \code{ggplot} objects.
#'
#' @description
#' Generates smooth occupancy response curves with credible intervals for
#' each predictor. Supports optional faceting by year and class.
#'
#' @author Thiago Bicudo
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line labs theme_minimal theme
#' @importFrom ggplot2 element_text element_rect
#' @importFrom dplyr group_split group_keys
#'
#' @export
plot_occ_predictions <- function(
    pred_occ_df,
    type_levels = NULL,
    color_ribbon = "#CCCCCC",
    color_line = "#1B3A5C",
    color_text = "#1B3A5C",
    color_strip_bg = "#E8EEF4",
    ribbon_alpha = 0.55,
    strip_class_prefix = ""
) {

  # ── Input validation ───────────────────────────────────────────────
  stopifnot(
    "pred_occ_df must be a data.frame" = is.data.frame(pred_occ_df),
    "Required columns missing" =
      all(c("variable", "value", "mean", "criL", "criU") %in% names(pred_occ_df))
  )

  # ── Optional cleaning of class labels ──────────────────────────────
  if ("class" %in% names(pred_occ_df) && strip_class_prefix != "") {
    pred_occ_df$class <- gsub(strip_class_prefix, "", pred_occ_df$class)
  }

  # ── Define grouping structure ──────────────────────────────────────
  grouping_vars <- c("variable")

  if ("year" %in% names(pred_occ_df)) {
    grouping_vars <- c(grouping_vars, "year")
  }

  if (!is.null(type_levels) && "class" %in% names(pred_occ_df)) {
    grouping_vars <- c(grouping_vars, "class")
  }

  # Split data
  split_data <- pred_occ_df %>%
    dplyr::group_by(across(all_of(grouping_vars))) %>%
    dplyr::group_split()

  keys <- pred_occ_df %>%
    dplyr::group_by(across(all_of(grouping_vars))) %>%
    dplyr::group_keys()

  # ── Plot builder ───────────────────────────────────────────────────
  plot_list <- vector("list", length(split_data))

  for (i in seq_along(split_data)) {

    df <- split_data[[i]]
    key <- keys[i, , drop = FALSE]

    # Build title dynamically
    title_parts <- unlist(key)
    plot_title <- paste(title_parts, collapse = " | ")

    p <- ggplot2::ggplot(df, ggplot2::aes(x = value, y = mean)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = criL, ymax = criU),
        fill = color_ribbon,
        alpha = ribbon_alpha
      ) +
      ggplot2::geom_line(color = color_line, linewidth = 1) +
      ggplot2::labs(
        x = unique(df$variable),
        y = "Occupancy probability",
        title = plot_title
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        text = ggplot2::element_text(color = color_text),
        plot.title = ggplot2::element_text(face = "bold"),
        strip.background = ggplot2::element_rect(fill = color_strip_bg)
      )

    plot_list[[i]] <- p
  }

  # ── Naming plots ───────────────────────────────────────────────────
  names(plot_list) <- apply(keys, 1, paste, collapse = "_")

  return(plot_list)
}
