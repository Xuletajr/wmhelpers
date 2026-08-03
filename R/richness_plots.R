# =============================================================================
# Richness visualization functions
# File: R/richness_plots.R
# =============================================================================


#' Categorical dot plot of species richness
#'
#' Generates a jittered dot plot of site-level richness estimates grouped by
#' a categorical variable, with a summary \code{pointrange} overlay showing
#' the mean ± SD per group. Suitable for comparing richness across regions,
#' habitat types, or other discrete groupings.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @param richness_data A \code{tibble} output from
#'   \code{build_richness_summary()}, containing at minimum columns \code{mean}
#'   and the variable named in \code{x_var}.
#' @param x_var Character. Column name for the x-axis (categorical grouping
#'   variable, e.g., \code{"region"}).
#' @param y_var Character. Column name for the y-axis. Defaults to
#'   \code{"mean"}.
#' @param group_var Optional character. Column for optional dodging and legend.
#'   Defaults to \code{NULL}.
#' @param color_var Optional character. Column name for color/fill mapping.
#'   If \code{NULL}, uses \code{color_point} for all points. Defaults to
#'   \code{NULL}.
#' @param title Character or \code{NULL}. Plot title. Defaults to \code{NULL}.
#' @param x_label Character or \code{NULL}. X-axis label. If \code{NULL},
#'   uses \code{x_var}. Defaults to \code{NULL}.
#' @param y_label Character. Y-axis label. Defaults to
#'   \code{"Estimated Species Richness"}.
#' @param dodge_width Numeric. Width for position dodging. Defaults to
#'   \code{0.2}.
#' @param color_point Character. Hex color for individual points when
#'   \code{color_var} is \code{NULL}. Defaults to \code{"#008AA1"}.
#' @param color_summary Character. Hex color for the summary pointrange
#'   overlay. Defaults to \code{"#1B3644"}.
#' @param color_fill Character. Hex fill color for the summary diamond.
#'   Defaults to \code{"#1B3644"}.
#' @param text_color Character. Hex color for axis text. Defaults to
#'   \code{"#1B3644"}.
#'
#' @return A \code{ggplot2} object.
#'
#' @seealso \code{\link{build_richness_summary}}
#'
#' @importFrom ggplot2 ggplot aes geom_point stat_summary labs theme_bw theme element_text element_blank scale_color_manual position_jitter
#'
#' @export
plot_richness_categorical <- function(richness_data,
                                      x_var,
                                      y_var       = "mean",
                                      group_var   = NULL,
                                      color_var   = NULL,
                                      title       = NULL,
                                      x_label     = NULL,
                                      y_label     = "Estimated Species Richness",
                                      dodge_width = 0.2,
                                      color_point   = "#008AA1",
                                      color_summary = "#1B3644",
                                      color_fill    = "#1B3644",
                                      text_color    = "#1B3644") {

  stopifnot(
    "richness_data must be a data frame" = is.data.frame(richness_data),
    "x_var must be a single string"      = is.character(x_var) && length(x_var) == 1,
    "x_var must exist in richness_data"  = x_var %in% names(richness_data),
    "y_var must exist in richness_data"  = y_var %in% names(richness_data)
  )

  p <- ggplot2::ggplot(richness_data,
                       ggplot2::aes(x = .data[[x_var]],
                                    y = .data[[y_var]]))

  if (!is.null(color_var)) {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(color = .data[[color_var]]),
        size     = 3,
        position = ggplot2::position_jitter(width = 0.1, height = 0),
        alpha    = 0.8
      )
  } else {
    p <- p +
      ggplot2::geom_point(
        size     = 3,
        color    = color_point,
        position = ggplot2::position_jitter(width = 0.1, height = 0),
        alpha    = 0.8
      )
  }

  p <- p +
    ggplot2::stat_summary(
      geom     = "pointrange",
      fun.data = "mean_sdl",
      fun.args = list(mult = 1),
      size     = 1,
      color    = color_summary,
      shape    = 23,
      fill     = color_fill
    ) +
    ggplot2::labs(
      title = title,
      x     = x_label %||% x_var,
      y     = y_label
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title    = ggplot2::element_text(size = 18, color = text_color),
      axis.text     = ggplot2::element_text(size = 14, color = text_color),
      panel.grid    = ggplot2::element_blank(),
      legend.position = "bottom"
    )

  if (!is.null(color_var)) {
    p <- p + ggplot2::scale_color_viridis_d(option = "mako",
                                            begin = 0.15, end = 0.9)
  }

  p
}


#' Ridgeline density plot of species richness by group
#'
#' Generates a ridgeline (joyplot) density plot of posterior species richness
#' samples grouped by a categorical variable (e.g., region, landcover,
#' habitat). Y-axis position of each group is ordered by median richness
#' (low to high), producing a visual gradient. Fill color is assigned
#' independently of that ordering, so that a given category keeps the same
#' color across different plots or data subsets.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @param richness_density A long-format data frame, typically the output of
#'   \code{build_richness_density()}. Must contain the columns \code{richness}
#'   and the column named in \code{group_var}.
#' @param group_var Character. Name of the grouping column in
#'   \code{richness_density} (e.g., \code{"landcover"}, \code{"region"}).
#' @param show_median Logical. If \code{TRUE} (default), draws a vertical
#'   line at the median of each group's density.
#' @param alpha Numeric. Fill transparency for the density ridges.
#'   Defaults to \code{0.7}.
#' @param scale Numeric. Overlap scale between ridges, passed to
#'   \code{ggridges::geom_density_ridges()}. Defaults to \code{1.5}.
#' @param palette Character vector of 2 or more hex colors used as anchor
#'   points for interpolation. Defaults to
#'   \code{c("#008AA1", "#B90078", "#00D973")}. Ignored if \code{color_map}
#'   is supplied.
#' @param color_map Optional named character vector mapping specific category
#'   names (matching the values in \code{group_var}) to specific hex colors.
#'   Use this for full manual control (e.g., a fixed color convention across
#'   projects). If supplied, takes priority over \code{palette} and every
#'   category present in the data must have an entry in \code{color_map}.
#'   Defaults to \code{NULL}.
#' @param title Character or \code{NULL}. Plot title. Defaults to \code{NULL}.
#' @param x_label Character. X-axis label. Defaults to
#'   \code{"Estimated Species Richness"}.
#' @param y_label Character or \code{NULL}. Y-axis label. Defaults to
#'   \code{NULL}.
#' @param text_color Character. Hex color for axis text and titles.
#'   Defaults to \code{"#1B3644"}.
#'
#' @return A \code{ggplot2} object.
#'
#' @details
#' Two independent orderings are used:
#' \itemize{
#'   \item \strong{Y-axis position}: groups are ordered by median richness
#'     (ascending), producing the classic ridgeline visual gradient from
#'     low to high richness.
#'   \item \strong{Fill color}: when \code{color_map} is \code{NULL},
#'     colors are interpolated from \code{palette} and assigned by
#'     \strong{alphabetical order} of category names, not by median rank.
#'     Because a category's name doesn't change between datasets or
#'     subsets, its color stays the same across different plots even when
#'     its median richness rank shifts — solving the common problem of the
#'     same category (e.g., \code{"Grassland"}) receiving different colors
#'     in different figures.
#' }
#'
#' Use \code{color_map} instead of \code{palette} only when you need to
#' override specific categories with fixed brand or convention colors
#' (e.g., always green for forest, always blue for water) regardless of
#' alphabetical order. The two arguments are not meant to be combined:
#' if \code{color_map} is supplied, \code{palette} is ignored entirely.
#'
#' @seealso \code{\link{build_richness_density}}, \code{\link{extract_richness}}
#'
#' @importFrom ggplot2 ggplot aes labs theme_bw theme element_text element_line
#' @importFrom ggplot2 element_blank scale_fill_manual after_scale
#' @importFrom ggridges geom_density_ridges
#' @importFrom colorspace darken
#' @importFrom grDevices colorRampPalette
#'
#' @export
plot_richness_ridgeline <- function(richness_density,
                                    group_var,
                                    show_median = TRUE,
                                    alpha       = 0.7,
                                    scale       = 1.5,
                                    palette     = c("#008AA1", "#B90078", "#00D973"),
                                    color_map   = NULL,
                                    title       = NULL,
                                    x_label     = "Estimated Species Richness",
                                    y_label     = NULL,
                                    text_color  = "#1B3644") {

  # ── Input validation ─────────────────────────────────────────────────────────
  stopifnot(
    "richness_density must be a data frame" = is.data.frame(richness_density),
    "group_var must be a single string"     = is.character(group_var) &&
      length(group_var) == 1,
    "group_var must exist in richness_density" = group_var %in%
      names(richness_density),
    "'richness' column must exist"          = "richness" %in%
      names(richness_density)
  )

  # ── Y-axis order: by median richness (ascending) ─────────────────────────────
  group_order <- richness_density %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::summarise(median_r = stats::median(richness), .groups = "drop") %>%
    dplyr::arrange(median_r) %>%
    dplyr::pull(.data[[group_var]])

  richness_density <- richness_density %>%
    dplyr::mutate(!!group_var := factor(.data[[group_var]], levels = group_order))

  n_groups <- length(group_order)

  # ── Color assignment ──────────────────────────────────────────────────────────
  if (!is.null(color_map)) {

    # ── Manual mode: explicit name -> color dictionary ─────────────────────────
    missing_colors <- setdiff(group_order, names(color_map))
    if (length(missing_colors) > 0) {
      stop(
        "color_map is missing colors for: ",
        paste(missing_colors, collapse = ", ")
      )
    }

    palette_final <- color_map[group_order]

  } else {

    # ── Automatic mode: interpolate from palette, assign alphabetically ────────
    # Colors are tied to category NAME (alphabetical order), not to median
    # rank, so a given category keeps the same color across different plots
    # or data subsets, regardless of where it falls in the richness gradient.
    alpha_order   <- sort(group_order)
    alpha_colors  <- (grDevices::colorRampPalette(palette))(n_groups)
    color_by_name <- stats::setNames(alpha_colors, alpha_order)

    palette_final <- color_by_name[group_order]
  }

  # ── Plot ──────────────────────────────────────────────────────────────────────
  ggplot2::ggplot(richness_density,
                  ggplot2::aes(x = richness, y = .data[[group_var]],
                               fill = .data[[group_var]])) +
    ggridges::geom_density_ridges(
      ggplot2::aes(color = ggplot2::after_scale(colorspace::darken(fill, 0.4))),
      alpha = alpha, scale = scale, linewidth = 1,
      quantile_lines = show_median, quantiles = 0.5
    ) +
    ggridges::geom_density_ridges(
      ggplot2::aes(vline_color = ggplot2::after_scale(colorspace::darken(fill, 0.5))),
      alpha = 0, scale = scale, linewidth = 0,
      vline_width = 1.2, quantile_lines = show_median, quantiles = 0.5
    ) +
    ggplot2::scale_fill_manual(values = palette_final) +
    ggplot2::labs(title = title, x = x_label, y = y_label) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title         = ggplot2::element_text(size = 18, color = text_color),
      axis.text          = ggplot2::element_text(size = 14, color = text_color),
      panel.grid.major.x = ggplot2::element_line(color = "grey90"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.position    = "none"
    )
}


#' Scatter plot of richness against a continuous covariate
#'
#' Generates a scatter plot of site-level mean richness estimates against a
#' continuous environmental covariate, with an optional smoothing trend line.
#' Suitable for exploring richness gradients along environmental axes.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @param richness_data A \code{tibble} output from
#'   \code{build_richness_summary()}, containing at minimum a \code{mean}
#'   column and the variable named in \code{x_var}.
#' @param x_var Character. Column name for the continuous covariate (x-axis).
#' @param color_var Optional character. Column name for color mapping (e.g.,
#'   \code{"region"}). If \code{NULL}, all points use \code{color_point}.
#'   Defaults to \code{NULL}.
#' @param smooth_method Character. Smoothing method passed to
#'   \code{geom_smooth()}. One of \code{"lm"}, \code{"loess"}, or \code{"gam"}.
#'   Defaults to \code{"lm"}.
#' @param ci Logical. If \code{TRUE}, shows a confidence interval around the
#'   trend line. Defaults to \code{TRUE}.
#' @param palette Optional character vector of colors for manual color scale
#'   when \code{color_var} is provided. If \code{NULL} or insufficient, falls
#'   back to the \code{viridis} \code{"mako"} palette. Defaults to \code{NULL}.
#' @param title Character or \code{NULL}. Plot title. Defaults to \code{NULL}.
#' @param x_label Character or \code{NULL}. X-axis label. If \code{NULL},
#'   uses \code{x_var}. Defaults to \code{NULL}.
#' @param y_label Character. Y-axis label. Defaults to
#'   \code{"Estimated Species Richness (mean)"}.
#' @param color_point Character. Hex color for points when \code{color_var}
#'   is \code{NULL}. Defaults to \code{"#008AA1"}.
#' @param color_smooth Character. Hex color for the smooth line and CI fill.
#'   Defaults to \code{"#1B3644"}.
#' @param text_color Character. Hex color for axis text. Defaults to
#'   \code{"#1B3644"}.
#'
#' @return A \code{ggplot2} object.
#'
#' @seealso \code{\link{build_richness_summary}}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_bw theme element_text element_blank scale_color_manual scale_color_viridis_d
#'
#' @export
plot_richness_continuous <- function(richness_data,
                                     x_var,
                                     color_var     = NULL,
                                     smooth_method = "lm",
                                     ci            = TRUE,
                                     palette       = NULL,
                                     title         = NULL,
                                     x_label       = NULL,
                                     y_label       = "Estimated Species Richness (mean)",
                                     color_point   = "#008AA1",
                                     color_smooth  = "#1B3644",
                                     text_color    = "#1B3644") {

  stopifnot(
    "richness_data must be a data frame" = is.data.frame(richness_data),
    "x_var must be a single string"      = is.character(x_var) && length(x_var) == 1,
    "x_var must exist in richness_data"  = x_var %in% names(richness_data),
    "'mean' column must exist"           = "mean" %in% names(richness_data)
  )

  p <- ggplot2::ggplot(richness_data,
                       ggplot2::aes(x = .data[[x_var]], y = mean))

  if (!is.null(color_var)) {

    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(color = .data[[color_var]]),
        size = 2.5, alpha = 0.8
      )

    n_groups <- length(unique(richness_data[[color_var]]))

    if (!is.null(palette) && length(palette) >= n_groups) {
      p <- p + ggplot2::scale_color_manual(values = palette)
    } else {
      p <- p + ggplot2::scale_color_viridis_d(option = "mako",
                                              begin = 0.15, end = 0.9)
    }

  } else {
    p <- p +
      ggplot2::geom_point(color = color_point, size = 2.5, alpha = 0.8)
  }

  p +
    ggplot2::geom_smooth(
      method    = smooth_method,
      se        = ci,
      color     = color_smooth,
      fill      = color_smooth,
      alpha     = 0.15,
      linewidth = 1.0
    ) +
    ggplot2::labs(
      title = title,
      x     = x_label %||% x_var,
      y     = y_label,
      color = color_var
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title      = ggplot2::element_text(size = 16, color = text_color),
      axis.text       = ggplot2::element_text(size = 13, color = text_color),
      panel.grid      = ggplot2::element_blank(),
      legend.position = if (!is.null(color_var)) "bottom" else "none"
    )
}
