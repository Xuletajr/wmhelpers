#' Build Occupancy Prediction Plot List
#'
#' Produces a named list of \code{ggplot2} objects showing occupancy
#' probability responses across predictor gradients. Each plot corresponds
#' to a combination of variable \eqn{\times} year \eqn{\times} (optionally) class.
#'
#' @param pred_occ_df Tidy data frame returned by
#'   \code{\link{build_occ_summary}}. Must include columns:
#'   \code{variable}, \code{expanded_variable}, \code{seq_data},
#'   \code{mean}, \code{ciL}, \code{ciU}, and optionally
#'   \code{year} and \code{class}.
#' @param type_levels Optional character vector of class/type levels
#'   (e.g. habitat types). Defaults to \code{NULL}.
#' @param y_label Character. Label for the y-axis. Defaults to
#'   \code{"Occupancy probability"}. Override with
#'   \code{"Detection probability"} when plotting output from
#'   \code{\link{compute_det_predictions}}.
#' @param color_ribbon Character. Ribbon fill colour.
#'   Defaults to \code{"#CCCCCC"}.
#' @param color_line Character. Line colour.
#'   Defaults to \code{"#1B3A5C"}.
#' @param color_text Character. Text/axis colour.
#'   Defaults to \code{"#1B3A5C"}.
#' @param color_strip_bg Character. Facet strip background colour.
#'   Defaults to \code{"#E8EEF4"}.
#' @param ribbon_alpha Numeric. Transparency of ribbon.
#'   Defaults to \code{0.55}.
#' @param x_size Numeric. Font size for x-axis text and title.
#'   Defaults to \code{12}.
#' @param y_size Numeric. Font size for y-axis text and title.
#'   Defaults to \code{12}.
#' @param title_size Numeric. Font size for the plot title.
#'   Defaults to \code{13}.
#' @param legend_title_size Numeric. Font size for legend title text.
#'   Defaults to \code{12}.
#' @param legend_text_size Numeric. Font size for legend key labels.
#'   Defaults to \code{11}.
#' @param strip_class_prefix Character string to remove from class labels.
#'   Defaults to \code{""}.
#'
#' @return Named list of \code{ggplot} objects. Plot names are built from
#'   non-NA grouping keys (e.g. \code{"scosa"} for single-period models,
#'   \code{"scosa_2023"} for multi-year models).
#'
#' @description
#' Generates smooth response curves with credible intervals for
#' each predictor. Supports optional faceting by year and class.
#' Column names are aligned with the output of
#' \code{\link{build_occ_summary}}: x-axis uses \code{seq_data}
#' and credible interval bounds use \code{ciL}/\code{ciU}.
#'
#' When \code{year} or \code{class} are \code{NA} (e.g. single-period
#' models), they are automatically suppressed from plot titles and filenames.
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
    type_levels        = NULL,
    y_label             = "Occupancy probability",
    color_ribbon        = "#CCCCCC",
    color_line          = "#1B3A5C",
    color_text          = "#1B3A5C",
    color_strip_bg      = "#E8EEF4",
    ribbon_alpha        = 0.55,
    x_size              = 12,
    y_size              = 12,
    title_size          = 13,
    legend_title_size   = 12,
    legend_text_size    = 11,
    strip_class_prefix  = ""
) {

  # Input validation
  stopifnot(
    "pred_occ_df must be a data.frame" = is.data.frame(pred_occ_df),
    "Required columns missing" =
      all(c("variable", "expanded_variable", "seq_data",
            "mean", "ciL", "ciU") %in% names(pred_occ_df))
  )

  # Optional cleaning of class labels
  if ("class" %in% names(pred_occ_df) && strip_class_prefix != "") {
    pred_occ_df$class <- gsub(strip_class_prefix, "", pred_occ_df$class)
  }

  # Define grouping structure
  grouping_vars <- c("variable")

  if ("year" %in% names(pred_occ_df)) {
    grouping_vars <- c(grouping_vars, "year")
  }

  if (!is.null(type_levels) && "class" %in% names(pred_occ_df)) {
    grouping_vars <- c(grouping_vars, "class")
  }

  # Split data
  split_data <- pred_occ_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>%
    dplyr::group_split()

  keys <- pred_occ_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>%
    dplyr::group_keys()

  # Helper: clean key parts (remove NA and "variable" column)
  clean_key_parts <- function(key_row, exclude_col = "variable") {
    parts <- unlist(key_row)
    parts <- parts[!is.na(parts) & parts != "NA"]
    parts <- parts[names(parts) != exclude_col]
    parts
  }

  # Plot builder
  plot_list <- vector("list", length(split_data))

  for (i in seq_along(split_data)) {

    df  <- split_data[[i]]
    key <- keys[i, , drop = FALSE]

    # Title: non-NA parts excluding variable name
    title_parts <- clean_key_parts(key, exclude_col = "variable")
    plot_title  <- if (length(title_parts) == 0) "" else paste(title_parts, collapse = " | ")

    # x-axis label: use expanded_variable (human-readable) if available
    x_label <- unique(df$expanded_variable)
    if (length(x_label) != 1) x_label <- unique(df$variable)[1]

    p <- ggplot2::ggplot(df, ggplot2::aes(x = seq_data, y = mean)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = ciL, ymax = ciU),
        fill  = color_ribbon,
        alpha = ribbon_alpha
      ) +
      ggplot2::geom_line(color = color_line, linewidth = 1) +
      ggplot2::labs(
        x     = x_label,
        y     = y_label,
        title = plot_title
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        text              = ggplot2::element_text(color = color_text),
        plot.title        = ggplot2::element_text(face = "bold",
                                                  size = title_size),
        strip.background  = ggplot2::element_rect(fill = color_strip_bg),
        axis.title.x      = ggplot2::element_text(size = x_size),
        axis.text.x       = ggplot2::element_text(size = x_size),
        axis.title.y      = ggplot2::element_text(size = y_size),
        axis.text.y       = ggplot2::element_text(size = y_size),
        legend.title      = ggplot2::element_text(size = legend_title_size),
        legend.text       = ggplot2::element_text(size = legend_text_size)
      )

    plot_list[[i]] <- p
  }

  # Naming plots: use all non-NA key parts (including variable)
  names(plot_list) <- apply(keys, 1, function(x) {
    parts <- x[!is.na(x) & x != "NA"]
    paste(parts, collapse = "_")
  })

  return(plot_list)
}
