# =============================================================================
# Predicted occupancy map function
# File: R/plot_occupancy_map.R
# =============================================================================


#' Plot Predicted Occupancy Map for a Single Species
#'
#' Generates and saves a spatial map of predicted occupancy probability
#' (median posterior) for one species across an H3 hexagonal grid. Detection
#' points are overlaid as presence markers. Supports both batch production
#' runs (via a species index loop) and interactive inspection.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta") and Thiago Bicudo
#'
#' @param sp Integer. Index of the species to plot, relative to
#'   \code{species_name}.
#' @param psi_pred_quantiles A 3-dimensional array of predicted occupancy
#'   quantiles with dimensions \code{[quantile x species x sites]}.
#'   Dimension names must include \code{"median"} in the first dimension and
#'   species names in the second. Typically the output of
#'   \code{spOccupancy::predict.intMsPGOcc()} after quantile extraction.
#' @param hexdat An \code{sf} object representing the H3 hexagonal prediction
#'   grid. Must contain a column \code{h3_address} for joining with
#'   \code{psi_pred_quantiles}.
#' @param points An \code{sf} object or data frame of detection points. Must
#'   contain columns \code{lon} and \code{lat} and one column per species
#'   named with the scientific name (1 = detected, 0 = not detected).
#' @param species_name Character vector of scientific names, aligned with the
#'   second dimension of \code{psi_pred_quantiles}.
#' @param common_name Character vector of common names, aligned with
#'   \code{species_name}.
#' @param output_dir Character. Directory where PNG files are saved.
#' @param project_name Character. String appended to output filenames.
#' @param quantile_stat Character. Which quantile to map. Must be one of the
#'   values in \code{dimnames(psi_pred_quantiles)[[1]]}. Defaults to
#'   \code{"median"}.
#' @param color_gradient Character vector of 3 colors used for the fill
#'   gradient (low → mid → high occupancy probability). Defaults to
#'   \code{c("#F7F7F7", "#4393C3", "#053061")}.
#' @param color_point_fill Character. Fill color for presence points.
#'   Defaults to \code{"#D95F02"}.
#' @param color_point_stroke Character. Stroke color for presence points.
#'   Defaults to \code{"#1B3644"}.
#' @param point_size Numeric. Size of presence point markers.
#'   Defaults to \code{2.3}.
#' @param point_alpha Numeric. Transparency of presence point markers.
#'   Defaults to \code{0.5}.
#' @param width Numeric. Output image width in cm. Defaults to \code{35}.
#' @param height Numeric. Output image height in cm. Defaults to \code{30}.
#' @param dpi Numeric. Output resolution. Defaults to \code{600}.
#' @param print_plot Logical. If \code{TRUE}, prints the plot to the active
#'   graphics device. Defaults to \code{FALSE}.
#'
#' @return Invisibly returns the \code{ggplot2} object. The plot is also
#'   saved as a PNG file at
#'   \code{output_dir/<scientific_name_underscored>_<project_name>.png}.
#'
#' @details
#' The function pre-joins the predicted occupancy quantile for the requested
#' species to \code{hexdat} via \code{h3_address}, then overlays presence
#' points filtered from \code{points} where the species column equals 1.
#'
#' For batch production runs, call within a loop over species indices:
#' \preformatted{
#' for (sp in seq_along(species_name)) {
#'   plot_predictions(
#'     sp                = sp,
#'     psi_pred_quantiles = psi.pred.quantiles,
#'     hexdat            = hex_sf,
#'     points            = site_pts,
#'     species_name      = species_name,
#'     common_name       = common_name,
#'     output_dir        = paste0(config$output_maps, "MSOM/Species_audio/"),
#'     project_name      = config$project_name
#'   )
#' }
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_sf geom_point scale_fill_gradientn
#' @importFrom ggplot2 scale_shape_manual guide_legend guides guide_colorbar
#' @importFrom ggplot2 labs theme_bw theme element_blank element_text unit
#' @importFrom ggplot2 ggsave
#' @importFrom dplyr left_join filter select mutate
#' @importFrom tibble tibble
#' @importFrom ggspatial annotation_north_arrow annotation_scale north_arrow_fancy_orienteering
#'
#' @export
plot_occupancy_map <- function(sp,
                             psi_pred_quantiles,
                             hexdat,
                             points,
                             species_name,
                             common_name,
                             output_dir,
                             project_name,
                             quantile_stat      = "median",
                             color_gradient     = c("#F7F7F7", "#4393C3", "#053061"),
                             color_point_fill   = "#D95F02",
                             color_point_stroke = "#1B3644",
                             point_size         = 2.3,
                             point_alpha        = 0.5,
                             width              = 35,
                             height             = 30,
                             dpi                = 600,
                             print_plot         = FALSE) {

  # ── Input validation ─────────────────────────────────────────────────────────
  stopifnot(
    "sp must be a single positive integer"          = is.numeric(sp) &&
      length(sp) == 1 &&
      sp >= 1,
    "psi_pred_quantiles must be a 3D array"         = is.array(psi_pred_quantiles) &&
      length(dim(psi_pred_quantiles)) == 3,
    "hexdat must be an sf object"                   = inherits(hexdat, "sf"),
    "hexdat must contain 'h3_address' column"       = "h3_address" %in% names(hexdat),
    "species_name must be a character vector"       = is.character(species_name),
    "common_name must be a character vector"        = is.character(common_name),
    "species_name and common_name must have same length" =
      length(species_name) == length(common_name),
    "sp index out of range"                         = sp <= length(species_name),
    "quantile_stat must be in dimnames[[1]]"        = quantile_stat %in%
      dimnames(psi_pred_quantiles)[[1]],
    "output_dir must be a single string"            = is.character(output_dir) &&
      length(output_dir) == 1,
    "project_name must be a single string"          = is.character(project_name) &&
      length(project_name) == 1
  )

  # Normalise output_dir
  if (!endsWith(output_dir, "/")) output_dir <- paste0(output_dir, "/")

  # ── Identify species ─────────────────────────────────────────────────────────
  sp_name   <- species_name[sp]
  sp_common <- common_name[sp]

  if (!sp_name %in% dimnames(psi_pred_quantiles)[[2]]) {
    stop(
      "Species '", sp_name, "' not found in dimnames(psi_pred_quantiles)[[2]]."
    )
  }

  # ── Merge predicted occupancy to hex grid ────────────────────────────────────
  psi_tbl <- tibble::tibble(
    h3_address = dimnames(psi_pred_quantiles)[[3]],
    psi_stat   = psi_pred_quantiles[quantile_stat, sp_name, ]
  )

  psi_pred <- hexdat %>%
    dplyr::left_join(psi_tbl, by = "h3_address")

  # ── Filter presence points ───────────────────────────────────────────────────
  if (!sp_name %in% names(points)) {
    warning(
      "Column '", sp_name, "' not found in points. ",
      "No presence points will be plotted."
    )
    presence_pts <- data.frame(lon = numeric(0), lat = numeric(0))
  } else {
    presence_pts <- points %>%
      dplyr::filter(.data[[sp_name]] == 1) %>%
      dplyr::select(lon, lat)
  }

  # ── Build plot ───────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot() +

    # Hex grid fill
    ggplot2::geom_sf(
      data  = psi_pred,
      ggplot2::aes(fill = psi_stat),
      color = NA
    ) +

    # Presence points
    ggplot2::geom_point(
      data   = presence_pts,
      ggplot2::aes(x = lon, y = lat, shape = "Presence"),
      fill   = color_point_fill,
      colour = color_point_stroke,
      size   = point_size,
      stroke = 0.8,
      alpha  = point_alpha
    ) +

    # Fill color scale
    ggplot2::scale_fill_gradientn(
      colours = color_gradient,
      name    = "Probability of Occupancy",
      na.value = "grey92"
    ) +

    # Shape legend
    ggplot2::scale_shape_manual(
      name   = "Species detected",
      values = c("Presence" = 21),
      labels = c("Presence" = ""),
      guide  = ggplot2::guide_legend(
        override.aes = list(
          fill   = color_point_fill,
          colour = color_point_stroke,
          size   = 5,
          stroke = 0.9
        )
      )
    ) +

    # North arrow
    ggspatial::annotation_north_arrow(
      location    = "tl",
      which_north = "true",
      style       = ggspatial::north_arrow_fancy_orienteering()
    ) +

    # Scale bar
    ggspatial::annotation_scale(
      location = "bl",
      bar_cols = c("grey60", "white"),
      text_col = "black"
    ) +

    ggplot2::labs(
      title    = paste0("Predicted Occupancy \u2014 ", sp_common),
      subtitle = sp_name
    ) +

    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      axis.title       = ggplot2::element_blank(),
      panel.grid       = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(size = 16, face = "bold"),
      plot.subtitle    = ggplot2::element_text(face = "italic"),
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.box       = "horizontal",
      legend.key.size  = ggplot2::unit(2, "line"),
      legend.spacing.x = ggplot2::unit(0.3, "in"),
      legend.title     = ggplot2::element_text(size = 16)
    ) +

    # Legend order: presence points first, color bar second
    ggplot2::guides(
      fill  = ggplot2::guide_colorbar(order = 2),
      shape = ggplot2::guide_legend(order = 1)
    )

  # ── Print (optional) ─────────────────────────────────────────────────────────
  if (print_plot) print(p)

  # ── Save ─────────────────────────────────────────────────────────────────────
  filename <- paste0(
    output_dir,
    gsub(" ", "_", sp_name), "_", project_name, ".png"
  )

  message("Saving: ", sp_name)

  ggplot2::ggsave(
    filename = filename,
    plot     = p,
    width    = width,
    height   = height,
    dpi      = dpi,
    units    = "cm",
    bg       = "white"
  )

  invisible(p)
}
