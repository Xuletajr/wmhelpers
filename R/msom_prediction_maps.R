# =============================================================================
# File: msomHelper.R
#
# Description:
#   Helper functions for generating spatial prediction maps from a fitted
#   Hierarchical Multi-Species Occupancy Model (MSOM). Intended to be sourced
#   by project-specific scripts. All functions depend on a `config` object
#   loaded from a YAML project configuration file.
#
# Authors:
#   Original: Kris Harmon
#   Updated and expanded: Jose Ribeiro 'Xuleta' and Thiago Bicudo
#
# Date created: April, 2025
# Last updated: January, 2026
#
# Dependencies:
#   tidyverse, spOccupancy, sf, ggspatial, viridis, stars, ggtext, h3jsr, yaml
# =============================================================================


# -----------------------------------------------------------------------------
#' Load MSOM Input Data
#'
#' Loads the workspace RData file produced by the multispecies detection loop,
#' filters the species list to the target group, and returns cleaned vectors
#' of species codes, scientific names, and common names.
#'
#' @param config A named list loaded from the project YAML configuration file.
#'   Must contain elements `workspace_path` and `project_name`.
#' @param sp_codes Character vector of species codes to retain (e.g.,
#'   `sp_codes_frogs`).
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{spp_list_ea}{Filtered data frame with columns `species`, `group`,
#'       and `common_name`.}
#'     \item{species_name}{Character vector of scientific names.}
#'     \item{common_name}{Character vector of common names.}
#'   }
#'
#' @examples
#' \dontrun{
#' config    <- yaml::yaml.load_file("project_config.yaml")
#' msom_data <- load_msom_data(config, sp_codes_frogs)
#' }
#'
#' @export
load_msom_data <- function(config, sp_codes) {

  load(paste0(config$workspace_path,
              "data_for_MSOM_",
              config$project_name,
              ".RData"),
       envir = parent.frame())

  spp_list_ea <- get("spp_list_ea", envir = parent.frame())

  spp_list_ea <- spp_list_ea %>%
    #dplyr::filter(species %in% sp_codes) %>%
    dplyr::select(species, group, common_name)

  list(
    spp_list_ea  = spp_list_ea,
    species_name = spp_list_ea$species,
    common_name  = spp_list_ea$common_name
  )
}


# -----------------------------------------------------------------------------
#' Create Output Directories for MSOM Maps and Rasters
#'
#' Creates the required output directory structure for maps and rasters if
#' they do not already exist.
#'
#' @param config A named list loaded from the project YAML configuration file.
#'   Must contain elements `output_maps` and `output_rasters`.
#'
#' @return Invisibly returns the path to the base maps directory.
#'
#' @examples
#' \dontrun{
#' config <- yaml::yaml.load_file("project_config.yaml")
#' create_msom_dirs(config)
#' }
#'
#' @export
create_msom_dirs <- function(config) {

  base_dir <- paste0(config$output_maps, "Occupancy_Models/")

  dir.create(base_dir,
             recursive    = TRUE,
             showWarnings = TRUE)

  dir.create(paste0(config$output_rasters, "Occupancy_Models/Species/"),
             recursive    = TRUE,
             showWarnings = FALSE)

  invisible(base_dir)
}


# -----------------------------------------------------------------------------
#' Import MSOM Prediction Quantiles
#'
#' Reads the RDS files containing species-level occupancy probability quantiles
#' and community richness quantiles produced by the MSOM prediction chunks.
#'
#' @param config A named list loaded from the project YAML configuration file.
#'   Must contain elements `workspace_path` and `project_name`.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{psi.pred.quantiles}{3-D array \[quantile × species × site\] of
#'       occupancy probability quantiles.}
#'     \item{richness.quantiles}{Matrix \[site × quantile\] of richness
#'       quantiles.}
#'   }
#'
#' @examples
#' \dontrun{
#' config <- yaml::yaml.load_file("project_config.yaml")
#' preds  <- import_msom_predictions(config)
#' }
#'
#' @export
import_msom_predictions <- function(config) {

  psi.pred.quantiles <- readRDS(
    paste0(config$workspace_path,
           "psi_pred_quantiles_",
           config$project_name,
           ".rds")
  )

  richness.quantiles <- readRDS(
    paste0(config$workspace_path,
           "richness_pred_quantiles_",
           config$project_name,
           ".rds")
  )

  list(
    psi.pred.quantiles = psi.pred.quantiles,
    richness.quantiles = richness.quantiles
  )
}


# -----------------------------------------------------------------------------
#' Load Raw Detection Table
#'
#' Reads the raw species-by-site detection CSV produced by the occupancy models
#' pipeline, removes sites with no detections across all species, and subsets
#' to the modelled species codes.
#'
#' @param config A named list loaded from the project YAML configuration file.
#'   Must contain element `output_tables`.
#' @param sp_codes Character vector of species codes to retain.
#'
#' @return A data frame with columns `Name` (site identifier) followed by one
#'   column per species in `sp_codes`, coded as 0/1 detection.
#'
#' @examples
#' \dontrun{
#' config <- yaml::yaml.load_file("project_config.yaml")
#' allsp  <- load_detection_table(config, sp_codes_frogs)
#' }
#'
#' @export
load_detection_table <- function(config, sp_codes) {

  dir_path <- file.path(config$output_tables, "Occupancy_Models/")
  pattern  <- "raw_detections.*\\.csv"

  file_path <- file.path(
    dir_path,
    list.files(path = dir_path, pattern = pattern)
  )

  readr::read_csv(file_path) %>%
    dplyr::filter(!dplyr::if_all(.cols = 2:ncol(.), ~ .x == 0)) %>%
    dplyr::rename(Name = sites) %>%
    dplyr::select(Name, dplyr::all_of(sp_codes))
}


# -----------------------------------------------------------------------------
#' Load and Prepare Site Shapefile
#'
#' Reads the site-level point shapefile, standardises site names to lowercase,
#' optionally applies project-specific name corrections, and filters to sites
#' present in the detection table.
#'
#' @param config A named list loaded from the project YAML configuration file.
#'   Must contain element `layers`.
#' @param allsp Data frame returned by [load_detection_table()], used to filter
#'   sites.
#' @param name_corrections A named character vector where each **name** is the
#'   raw (lowercased) site identifier found in the shapefile and each **value**
#'   is the corrected identifier that matches the detection table. Pass `NULL`
#'   (the default) to skip this step entirely.
#'
#'   Example:
#'   ```r
#'   name_corrections <- c(
#'     "cu02_new"     = "cu02",
#'     "mo35_correct" = "mo35",
#'     "mo45_2"       = "mo45"
#'   )
#'   ```
#'
#' @return An `sf` object with columns `Name`, `latitude`, `longitude`, and the
#'   original geometry.
#'
#' @details
#'   Name standardisation always converts site names to lowercase before
#'   applying `name_corrections`. The corrections are resolved via a simple
#'   lookup: any name found in `names(name_corrections)` is replaced by its
#'   corresponding value; all other names are left unchanged. This is
#'   equivalent to a `dplyr::recode()` call and avoids project-specific logic
#'   inside the helper.
#'
#' @examples
#' \dontrun{
#' config <- yaml::yaml.load_file("project_config.yaml")
#' allsp  <- load_detection_table(config, sp_codes_frogs)
#'
#' # With name corrections
#' name_corrections <- c(
#'   "cu02_new"     = "cu02",
#'   "mo35_correct" = "mo35",
#'   "mo45_2"       = "mo45"
#' )
#' siteshp <- load_site_shapefile(config, allsp, name_corrections)
#'
#' # Without name corrections
#' siteshp <- load_site_shapefile(config, allsp)
#' }
#'
#' @export
load_site_shapefile <- function(config, allsp, name_corrections = NULL) {

  shp_path <- list.files(
    path       = file.path(config$layers, "Analysis_Layers"),
    pattern    = "^sitedata.*\\.shp$",
    full.names = TRUE
  )

  shp <- sf::st_read(shp_path) %>%
    dplyr::select(site_name, latitude, longitude) %>%
    dplyr::mutate(Name = stringr::str_to_lower(site_name))

  ## Apply project-specific name corrections if provided
  if (!is.null(name_corrections)) {
    if (!is.character(name_corrections) || is.null(names(name_corrections))) {
      stop(
        "`name_corrections` must be a named character vector, e.g.\n",
        '  c("old_name" = "new_name", "another_old" = "another_new")'
      )
    }
    shp <- shp %>%
      dplyr::mutate(
        Name = dplyr::recode(Name, !!!name_corrections)
      )
  }

  shp %>%
    dplyr::filter(Name %in% allsp$Name)
}


# -----------------------------------------------------------------------------
#' Build H3 Hexagonal Grid from Prediction Sites
#'
#' Converts the H3 cell addresses stored in the third dimension of the
#' occupancy prediction array into an `sf` polygon object.
#'
#' @param psi.pred.quantiles 3-D array \[quantile × species × site\] of
#'   occupancy probability quantiles. The third dimension names must be valid
#'   H3 cell addresses.
#'
#' @return An `sf` object with one polygon row per H3 hexagon and a column
#'   `h3_address` containing the cell identifiers.
#'
#' @examples
#' \dontrun{
#' hex_sf <- build_hex_grid(preds$psi.pred.quantiles)
#' }
#'
#' @export
build_hex_grid <- function(psi.pred.quantiles) {

  h3jsr::cell_to_polygon(
    dimnames(psi.pred.quantiles)[[3]],
    simple = FALSE
  )
}


# -----------------------------------------------------------------------------
#' Combine Site Coordinates with Detection Data
#'
#' Projects site locations to the CRS of the hexagonal grid, extracts UTM
#' coordinates, and joins with the detection table.
#'
#' @param siteshp `sf` object returned by [load_site_shapefile()].
#' @param hex_sf `sf` object returned by [build_hex_grid()]. Used only to
#'   obtain the target CRS.
#' @param allsp Data frame returned by [load_detection_table()].
#'
#' @return An `sf` object combining site names, UTM easting (`lon_utm`) and
#'   northing (`lat_utm`), per-species detection columns, and point geometry.
#'
#' @examples
#' \dontrun{
#' utm <- build_utm_points(siteshp, hex_sf, allsp)
#' }
#'
#' @export
build_utm_points <- function(siteshp, hex_sf, allsp) {

  siteshp %>%
    sf::st_transform(sf::st_crs(hex_sf)) %>%
    sf::st_centroid() %>%
    sf::st_coordinates() %>%
    data.frame() %>%
    dplyr::rename(lat_utm = Y, lon_utm = X) %>%
    cbind(siteshp, .) %>%
    dplyr::select(Name, lat_utm, lon_utm) %>%
    dplyr::left_join(allsp, by = "Name") %>%
    dplyr::select(1, 4:(ncol(.) - 1), 2:3, ncol(.))
}


# -----------------------------------------------------------------------------
#' Plot and Export Species Occupancy Prediction Map
#'
#' Generates a ggplot2 map showing the median posterior occupancy probability
#' on an H3 hexagonal grid, overlaid with detection points for a single
#' species, and saves it to disk as a PNG.
#'
#' @param sp Integer index into `species_name` and `common_name` vectors.
#' @param hexdat `sf` object returned by [build_hex_grid()].
#' @param points `sf` object returned by [build_utm_points()].
#' @param psi.pred.quantiles 3-D array \[quantile × species × site\].
#' @param species_name Character vector of scientific names (same order as
#'   model species dimension).
#' @param common_name Character vector of common names (same order).
#' @param config A named list loaded from the project YAML configuration file.
#'   Must contain `output_maps`, `project_name`, `color2`, and
#'   `color10`.
#'
#' @return Invisibly returns the `ggplot` object. Side effect: saves a PNG
#'   file to `<config$output_maps>/MSOM/Species/<sp_name>_<project>.png`.
#'
#' @examples
#' \dontrun{
#' plot_predictions(
#'   sp                = 1,
#'   hexdat            = hex_sf,
#'   points            = utm,
#'   psi.pred.quantiles = preds$psi.pred.quantiles,
#'   species_name      = msom_data$species_name,
#'   common_name       = msom_data$common_name,
#'   config            = config
#' )
#' }
#'
#' @export
plot_predictions <- function(sp, hexdat, points,
                             psi.pred.quantiles,
                             species_name, common_name,
                             config) {

  sp_name   <- species_name[sp]
  sp_common <- common_name[sp]

  ## Merge predicted occupancy onto hex grid
  psi_pred <- hexdat %>%
    dplyr::left_join(
      tibble::tibble(
        h3_address = dimnames(psi.pred.quantiles)$site,
        psi_median = psi.pred.quantiles["median", sp_name, ]
      ),
      by = "h3_address"
    )

  ## Filter to presence-only points for this species
  presence_pts <- points %>%
    dplyr::select(!!rlang::sym(sp_name), lon_utm, lat_utm) %>%
    dplyr::filter(!!rlang::sym(sp_name) == 1)

  ## Build plot
  p <- ggplot2::ggplot() +

    ggplot2::geom_sf(
      data  = psi_pred,
      ggplot2::aes(fill = psi_median),
      color = NA
    ) +

    ggplot2::geom_point(
      data = presence_pts,
      ggplot2::aes(x = lon_utm, y = lat_utm, shape = "Presence"),
      fill   = config$color10,
      colour = config$color2,
      size   = 2.3,
      stroke = 0.8,
      alpha  = 0.5
    ) +

    ggplot2::scale_fill_viridis_c(
      option    = "rocket",
      direction = -1,
      name      = "Probability of Occupancy"
    ) +

    ggplot2::scale_shape_manual(
      name   = "Species detected",
      values = c("Presence" = 21),
      labels = c("Presence" = ""),
      guide  = ggplot2::guide_legend(
        override.aes = list(
          fill   = config$color10,
          colour = config$color2,
          size   = 5,
          stroke = 0.9
        )
      )
    ) +

    ggspatial::annotation_north_arrow(
      location    = "tl",
      which_north = "true",
      style       = ggspatial::north_arrow_fancy_orienteering()
    ) +

    ggspatial::annotation_scale(
      location = "bl",
      bar_cols = c("grey60", "white"),
      text_col = "black"
    ) +

    ggplot2::labs(
      title    = paste0("Predicted Occupancy for ", sp_common),
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

    ggplot2::guides(
      fill  = ggplot2::guide_colorbar(order = 2),
      shape = ggplot2::guide_legend(order = 1,
                                    override.aes = list(size = 9))
    )

  invisible(p)

  ggplot2::ggsave(
    filename = paste0(
      config$output_maps,
      "Occupancy_Models/Species/",
      sp_name, "_", config$project_name, ".png"
    ),
    plot   = p,
    width  = 35,
    height = 30,
    dpi    = 600,
    units  = "cm",
    bg     = "white"
  )

  invisible(p)
}


# -----------------------------------------------------------------------------
#' Plot and Export Community Richness Map
#'
#' Generates a ggplot2 map of the posterior median community richness on the
#' H3 hexagonal grid and saves it to disk as a PNG.
#'
#' @param hex_sf `sf` object returned by [build_hex_grid()].
#' @param richness.quantiles Matrix \[site × quantile\] of richness quantiles.
#' @param config A named list loaded from the project YAML configuration file.
#'   Must contain `output_maps` and `project_name`.
#'
#' @return Invisibly returns the `ggplot` object. Side effect: saves a PNG
#'   file to
#'   `<config$output_maps>/Occupancy_Models/msom_median_richness_<project>.png`.
#'
#' @examples
#' \dontrun{
#' plot_richness_map(hex_sf, preds$richness.quantiles, config)
#' }
#'
#' @export
plot_richness_map <- function(hex_sf, richness.quantiles, config) {

  p <- hex_sf %>%
    dplyr::left_join(
      tibble::tibble(
        h3_address      = dimnames(richness.quantiles)$site,
        richness_median = richness.quantiles[, "median"]
      ),
      by = "h3_address"
    ) %>%

    ggplot2::ggplot() +

    ggplot2::geom_sf(
      ggplot2::aes(fill = richness_median),
      color = NA
    ) +

    ggplot2::scale_fill_viridis_c(
      option    = "rocket",
      direction = -1,
      name      = "Predicted Richness"
    ) +

    ggspatial::annotation_north_arrow(
      location    = "tl",
      which_north = "true",
      style       = ggspatial::north_arrow_fancy_orienteering()
    ) +

    ggspatial::annotation_scale(
      location = "bl",
      bar_cols = c("grey60", "white"),
      text_col = "black"
    ) +

    ggplot2::labs(
      title    = "Species richness",
      subtitle = "(posterior median)"
    ) +

    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      axis.title        = ggplot2::element_blank(),
      panel.grid        = ggplot2::element_blank(),
      axis.text         = ggplot2::element_blank(),
      axis.ticks        = ggplot2::element_blank(),
      plot.title        = ggplot2::element_text(size = 16, face = "bold"),
      plot.subtitle     = ggplot2::element_text(),
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.box        = "horizontal",
      legend.background = ggplot2::element_blank(),
      legend.text       = ggplot2::element_text(size = 14),
      legend.title      = ggplot2::element_text(size = 16, face = "bold"),
      legend.key.size   = ggplot2::unit(2, "line"),
      legend.spacing.x  = ggplot2::unit(0.3, "in")
    )

  print(p)

  ggplot2::ggsave(
    filename = normalizePath(
      paste0(config$output_maps,
             "Occupancy_Models/msom_median_richness_",
             config$project_name, ".png"),
      winslash = "/",
      mustWork = FALSE
    ),
    width  = 35,
    height = 30,
    dpi    = 600,
    units  = "cm",
    bg     = "white"
  )

  invisible(p)
}


# -----------------------------------------------------------------------------
#' Export Per-Species Occupancy Predictions as GeoTIFF Rasters
#'
#' Reprojects the H3 hexagonal grid to a metric CRS, creates a raster
#' template, rasterises the median occupancy prediction for each species, and
#' writes individual GeoTIFF files.
#'
#' @param hex_sf `sf` object returned by [build_hex_grid()].
#' @param psi.pred.quantiles 3-D array \[quantile × species × site\].
#' @param species_name Character vector of scientific names.
#' @param common_name Character vector of common names (used only for progress
#'   messages).
#' @param config A named list loaded from the project YAML configuration file.
#'   Must contain `output_rasters` and `project_name`.
#' @param crs_epsg Integer. EPSG code of the projected CRS to use for
#'   rasterisation. Default is `32161` (UTM Zone 19N / NAD83).
#' @param res_m Numeric. Raster resolution in metres. Default is `500`.
#'
#' @return Invisibly returns `NULL`. Side effect: writes one `.tif` file per
#'   species to `<config$output_rasters>/MSOM/`.
#'
#' @examples
#' \dontrun{
#' export_species_rasters(
#'   hex_sf             = hex_sf,
#'   psi.pred.quantiles = preds$psi.pred.quantiles,
#'   species_name       = msom_data$species_name,
#'   common_name        = msom_data$common_name,
#'   config             = config
#' )
#' }
#'
#' @export
export_species_rasters <- function(hex_sf,
                                   psi.pred.quantiles,
                                   species_name,
                                   common_name,
                                   config,
                                   crs_epsg = 32161,
                                   res_m    = 500) {

  hex_sf_utm <- sf::st_transform(hex_sf, crs_epsg)

  template <- stars::st_as_stars(
    sf::st_bbox(hex_sf_utm),
    dx     = res_m,
    dy     = res_m,
    values = NA
  )

  for (sp in seq_along(species_name)) {

    sp_name <- species_name[sp]

    psi_pred <- hex_sf_utm %>%
      dplyr::left_join(
        tibble::tibble(
          h3_address = dimnames(psi.pred.quantiles)$site,
          psi_median = psi.pred.quantiles["median", sp_name, ]
        ),
        by = "h3_address"
      )

    psi_raster <- stars::st_rasterize(
      psi_pred["psi_median"],
      template = template
    )

    out_path <- normalizePath(
      paste0(config$output_rasters,
             "Occupancy_Models/msom_",
             config$project_name,
             "_", sp_name, ".tif"),
      winslash = "/",
      mustWork = FALSE
    )

    stars::write_stars(psi_raster, out_path)
    message("Exported raster: ", sp_name)
  }

  invisible(NULL)
}

# END OF msomHelper.R ---------------------------------------------------------