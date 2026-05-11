# =============================================================================
# wmhelper — MSOM Spatial Prediction Functions
#
# Functions in this file (in pipeline order):
#
#   1. load_hexgrid_variables()     - read, clean and rename hexgrid CSV
#   2. prepare_hex_covariates()     - select model variables, scale, add intercept
#   3. project_coords()             - reproject lon/lat to a metric CRS
#   4. predict_occupancy_chunks()   - chunked predict() loop → psi + richness arrays
#   5. plot_species_occupancy_map() - quick single-species check map
#
# =============================================================================
# WHAT TO KEEP AS PLAIN SCRIPT (not wrapped in functions):
#   - setwd / rm / gc / library calls
#   - config loading, model and data loading
#   - out_bestmodel inspection calls (summary, $call, $x.p.names)
#   - chunk_size / n.omp.threads definitions  (project-specific tuning)
#   - saveRDS() calls
#   - commented-out code blocks
# =============================================================================


# =============================================================================
# 1. load_hexgrid_variables
# =============================================================================

#' Load and Clean Hexgrid Environmental Variables
#'
#' Reads the hexgrid covariate CSV from a directory, lowercases all column
#' names, removes rows with any \code{NA}, and applies a user-supplied
#' renaming map to align column names with those expected by the fitted model.
#'
#' @param env_path Character. Path to the directory containing the hexgrid
#'   CSV. The function auto-detects any file matching
#'   \code{hexdata*.csv} in that directory.
#' @param rename_map Optional named character vector where names are the
#'   \strong{new} column names and values are the \strong{original} column
#'   names (i.e. the same convention as \code{dplyr::rename}).
#'   Example: \code{c(footprint = "hfp2013pr", dist_water = "dist_water_log")}.
#'   Defaults to \code{NULL} (no renaming).
#' @param transform_map Optional named list of one-sided formulas or functions
#'   applied to specific columns after renaming. Names are column names,
#'   values are anonymous functions.
#'   Example: \code{list(dist_water = \(x) log1p(x))}.
#'   Defaults to \code{NULL}.
#'
#' @return A \code{tibble} with cleaned and renamed columns.
#'
#' @details
#' Auto-detection of the CSV file uses \code{list.files(pattern = "hexdata.*\\.csv")}.
#' If more than one file matches, the first match is used and a warning is
#' issued.
#'
#' @importFrom readr read_csv
#' @importFrom dplyr rename_with filter rename mutate
#' @importFrom stringr str_to_lower
#'
#' @export
load_hexgrid_variables <- function(env_path,
                                    rename_map     = NULL,
                                    transform_map  = NULL) {

  stopifnot(
    "env_path must be a character string" = is.character(env_path)
  )

  if (!endsWith(env_path, "/")) env_path <- paste0(env_path, "/")

  # ── Auto-detect CSV file ──────────────────────────────────────────────────────
  files <- list.files(env_path, pattern = "hexdata.*\\.csv", full.names = TRUE)

  if (length(files) == 0)
    stop("No file matching 'hexdata*.csv' found in: ", env_path, call. = FALSE)
  if (length(files) > 1)
    warning("Multiple hexdata files found — using: ", basename(files[1]),
            call. = FALSE)

  # ── Read and clean ─────────────────────────────────────────────────────────────
  dat <- readr::read_csv(files[1], show_col_types = FALSE) %>%
    dplyr::rename_with(stringr::str_to_lower) %>%
    dplyr::filter(!rowSums(is.na(.)))

  # ── Rename columns to match model names ───────────────────────────────────────
  if (!is.null(rename_map)) {
    stopifnot(
      "rename_map must be a named character vector" =
        is.character(rename_map) && !is.null(names(rename_map))
    )
    dat <- dplyr::rename(dat, !!!rename_map)
  }

  # ── Apply column transformations ──────────────────────────────────────────────
  if (!is.null(transform_map)) {
    for (col in names(transform_map)) {
      if (!col %in% names(dat))
        warning("transform_map column '", col, "' not found — skipping.",
                call. = FALSE)
      else
        dat[[col]] <- transform_map[[col]](dat[[col]])
    }
  }

  dat
}


# =============================================================================
# 2. prepare_hex_covariates
# =============================================================================

# =============================================================================
# 2.1. Extract Scaling Parameters from Training Sites
# =============================================================================

#' Extract Scaling Parameters from Training Sites
#'
#' Computes the mean and SD of each covariate from the training sites,
#' to be used later for scaling the hexagonal grid covariates.
#'
#' @param variables_of_interest Data frame of training-site covariates.
#' @param get_variables_of_interest Character vector of variable names to retain.
#' @param rename_vars Optional named character vector to rename columns before
#'   selecting. E.g. \code{c(agb = "agb_mean", slope = "slope_mean")}.
#'   Default is \code{NULL}.
#'
#' @return A data frame with columns \code{variable}, \code{mean}, and \code{sd}.
#'
#' @export
get_scaling_params <- function(variables_of_interest,
                               get_variables_of_interest,
                               rename_vars = NULL) {

  df <- variables_of_interest

  if (!is.null(rename_vars)) {
    df <- df %>% dplyr::rename(!!!rename_vars)
  }

  df %>%
    dplyr::select(dplyr::all_of(get_variables_of_interest)) %>%
    dplyr::summarise(dplyr::across(
      dplyr::everything(),
      list(
        mean = ~ mean(., na.rm = TRUE),
        sd   = ~ sd(.,   na.rm = TRUE)
      )
    )) %>%
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to      = c("variable", "stat"),
      names_pattern = "^(.+)_(mean|sd)$"
    ) %>%
    tidyr::pivot_wider(names_from = stat, values_from = value)
}

# =============================================================================
# 2.2 Prepare Hexagonal Grid Covariates for MSOM Prediction
# =============================================================================

#' Prepare Hexagonal Grid Covariates for MSOM Prediction
#'
#' Selects and scales the covariates needed by the best model from the
#' hexagonal grid data, using the training-site mean and SD to ensure
#' predictions are on the same scale as the fitted model.
#'
#' @param hexgrid_variables Data frame containing the hexagonal grid covariates.
#' @param out_bestmodel Fitted model object containing \code{x.names}.
#' @param scaling_params Data frame with columns \code{variable}, \code{mean},
#'   and \code{sd} — computed from the original training sites, not the hexgrid.
#' @param extra_columns Optional named list of additional columns to append to
#'   the design matrix (e.g. dummy or categorical variables). Default is \code{NULL}.
#' @param scale_vars Logical. If \code{TRUE} (default), scales covariates using
#'   \code{scaling_params}.
#'
#' @return A numeric matrix with an intercept column followed by the scaled
#'   model covariates, ready for prediction.
#'
#' @export
prepare_hex_covariates <- function(hexgrid_variables,
                                   out_bestmodel,
                                   scaling_params,
                                   extra_columns = NULL,
                                   scale_vars    = TRUE) {
  stopifnot(
    "hexgrid_variables must be a data frame"    = is.data.frame(hexgrid_variables),
    "out_bestmodel must have x.names"           = !is.null(out_bestmodel$x.names),
    "scaling_params must have variable/mean/sd" = all(c("variable", "mean", "sd") %in% names(scaling_params))
  )

  # Model covariates excluding the intercept
  model_vars   <- setdiff(out_bestmodel$x.names, "(Intercept)")
  missing_vars <- setdiff(model_vars, names(hexgrid_variables))

  if (length(missing_vars) > 0)
    stop("Columns missing from hexgrid_variables: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)

  # ── Select and optionally scale using training-site parameters ────────────
  hex_covs <- hexgrid_variables %>%
    dplyr::select(dplyr::all_of(model_vars))

  if (scale_vars) {
    hex_covs <- hex_covs %>%
      dplyr::mutate(dplyr::across(
        dplyr::where(is.numeric),
        ~ {
          var_name <- dplyr::cur_column()
          mu <- scaling_params$mean[scaling_params$variable == var_name]
          sg <- scaling_params$sd[scaling_params$variable == var_name]
          if (length(mu) == 0 || length(sg) == 0)
            stop("No scaling params found for variable: ", var_name, call. = FALSE)
          (. - mu) / sg
        }
      ))
  }

  # ── Sanity check ──────────────────────────────────────────────────────────
  message("Column means after scaling (should be ≈ 0):")
  print(round(colMeans(hex_covs), 5))
  message("Column SDs after scaling (should be ≈ 1):")
  print(round(apply(hex_covs, 2, sd), 5))

  # ── Build design matrix with intercept ───────────────────────────────────
  design <- data.frame("(Intercept)" = 1,
                       hex_covs,
                       check.names = FALSE)

  if (!is.null(extra_columns)) {
    for (col_name in names(extra_columns)) {
      design[[col_name]] <- extra_columns[[col_name]]
    }
  }

  as.matrix(design)
}


# =============================================================================
# 3. project_coords
# =============================================================================

#' Reproject Longitude/Latitude Coordinates to a Metric CRS
#'
#' Extracts longitude and latitude from a data frame, converts to an
#' \code{sf} point object in WGS84, reprojects to a user-specified metric
#' CRS, and returns a two-column tibble (\code{lon}, \code{lat}) of projected
#' coordinates.
#'
#' @param data Data frame containing coordinate columns.
#' @param lon_col Character. Name of the longitude column.
#'   Defaults to \code{"longitude"}.
#' @param lat_col Character. Name of the latitude column.
#'   Defaults to \code{"latitude"}.
#' @param target_crs Integer or character. EPSG code of the target projected
#'   CRS (e.g. \code{32161} for NAD83 Puerto Rico & Virgin Islands,
#'   \code{32618} for UTM zone 18N).
#'
#' @return A \code{tibble} with two numeric columns \code{lon} and \code{lat}
#'   in the projected coordinate system, one row per input point.
#'
#' @importFrom sf st_as_sf st_transform st_coordinates
#' @importFrom tibble as_tibble
#'
#' @export
project_coords <- function(data,
                            lon_col    = "longitude",
                            lat_col    = "latitude",
                            target_crs ) {

  stopifnot(
    "data must be a data frame"              = is.data.frame(data),
    "lon_col must be a column in data"       = lon_col %in% names(data),
    "lat_col must be a column in data"       = lat_col %in% names(data),
    "target_crs must be a single value"      = length(target_crs) == 1L
  )

  data %>%
    dplyr::select(dplyr::all_of(c(lon_col, lat_col))) %>%
    sf::st_as_sf(coords = c(lon_col, lat_col), crs = 4326) %>%
    sf::st_transform(target_crs) %>%
    sf::st_coordinates() %>%
    tibble::as_tibble() %>%
    setNames(c("lon", "lat"))
}


# =============================================================================
# 4. predict_occupancy_chunks
# =============================================================================

#' Predict Species Occupancy and Richness in Chunks
#'
#' Iterates over a covariate matrix and coordinate set in fixed-size chunks,
#' calling \code{predict()} on the fitted model for each chunk and storing
#' posterior quantiles of species occupancy (\code{psi}) and community
#' richness in pre-allocated arrays.
#'
#' Chunked prediction is necessary for large spatial grids that would exceed
#' available RAM if predicted in a single call.
#'
#' @param out_bestmodel Fitted \code{spOccupancy} model object.
#' @param hex_cov Numeric matrix as returned by
#'   \code{\link{prepare_hex_covariates}}.
#' @param coords Tibble or data frame with columns \code{lon} and \code{lat}
#'   as returned by \code{\link{project_coords}}.
#' @param site_ids Character vector of site identifiers (e.g. H3 addresses)
#'   used to label the site dimension of the output arrays.
#' @param sp_names Character vector of species names used to label the species
#'   dimension. Defaults to \code{out_bestmodel$sp.names}.
#' @param chunk_size Integer. Number of hexagons per prediction batch.
#'   Defaults to \code{1000}. Reduce if memory is limited.
#' @param n_omp_threads Integer. Number of OpenMP threads passed to
#'   \code{predict()}. Defaults to \code{1}.
#' @param ignore_re Logical. Passed as \code{ignore.RE} to \code{predict()}.
#'   Defaults to \code{FALSE}.
#' @param probs Numeric vector of length 3. Quantile probabilities for both
#'   psi and richness arrays. Defaults to \code{c(0.025, 0.5, 0.975)}.
#' @param quantile_names Character vector of length 3. Names for the quantile
#'   dimension. Defaults to \code{c("q025", "median", "q975")}.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\code{psi_quantiles}}{Three-dimensional array
#'     \code{[quantile × species × site]} of occupancy probability quantiles.}
#'   \item{\code{richness_quantiles}}{Matrix \code{[site × quantile]} of
#'     species richness quantiles.}
#' }
#'
#' @details
#' Progress is printed to the console as each chunk completes. Memory is
#' freed with \code{rm()} and \code{gc()} after each chunk.
#'
#' @export
predict_occupancy_chunks <- function(out_bestmodel,
                                      hex_cov,
                                      coords,
                                      site_ids,
                                      sp_names       = out_bestmodel$sp.names,
                                      chunk_size     = 1000L,
                                      n_omp_threads  = 1L,
                                      ignore_re      = FALSE,
                                      probs          = c(0.025, 0.5, 0.975),
                                      quantile_names = c("q025", "median", "q975")) {

  # ── Input validation ─────────────────────────────────────────────────────────
  stopifnot(
    "hex_cov must be a matrix"              = is.matrix(hex_cov),
    "coords must be a data frame"           = is.data.frame(coords),
    "nrow(hex_cov) must equal nrow(coords)" = nrow(hex_cov) == nrow(coords),
    "nrow(hex_cov) must equal length(site_ids)" =
      nrow(hex_cov) == length(site_ids),
    "probs must have length 3"              = length(probs) == 3,
    "quantile_names must have length 3"     = length(quantile_names) == 3
  )

  n_sites   <- nrow(hex_cov)
  n_species <- length(sp_names)
  n_chunks  <- ceiling(n_sites / chunk_size)

  # ── Pre-allocate output arrays ────────────────────────────────────────────────
  psi_quantiles      <- array(NA_real_,
                               dim = c(3L, n_species, n_sites))
  richness_quantiles <- array(NA_real_,
                               dim = c(n_sites, 3L))

  # ── Chunked prediction loop ───────────────────────────────────────────────────
  for (i in seq_len(n_chunks)) {

    cat("\n>>> Chunk", i, "of", n_chunks, "\n")

    start_row <- (i - 1L) * chunk_size + 1L
    end_row   <- min(i * chunk_size, n_sites)

    pred_chunk <- predict(
      out_bestmodel,
      hex_cov[start_row:end_row, , drop = FALSE],
      coords[start_row:end_row, ],
      ignore.RE      = ignore_re,
      n.omp.threads  = n_omp_threads,
      verbose        = FALSE
    )

    # Occupancy quantiles: quantile × species × site
    psi_quantiles[, , start_row:end_row] <- apply(
      pred_chunk$psi.0.samples, c(2, 3), quantile, probs = probs
    )

    ## ---- RICHNESS QUANTILES ---- ##
  # Sum across species (dim 2) for each MCMC x site -> [MCMC, sites]
  rich_mcmc <- apply(pred_chunk$z.0.samples, c(1, 3), sum)

  # Quantiles across MCMC (dim 1) per site -> [3, sites] -> t() -> [sites, 3]
  richness.quantiles[start_row:end_row, ] <- t(apply(
    rich_mcmc,
    2,
    quantile,
    probs = c(0.025, 0.5, 0.975)
  ))

    rm(pred_chunk, rich_mcmc)
    gc()
  }

  message("\nAll predictions complete.")

  # ── Assign dimension names ────────────────────────────────────────────────────
  dimnames(psi_quantiles) <- list(
    quantile = quantile_names,
    species  = sp_names,
    site     = site_ids
  )

  dimnames(richness_quantiles) <- list(
    site     = site_ids,
    quantile = quantile_names
  )

  list(
    psi_quantiles      = psi_quantiles,
    richness_quantiles = richness_quantiles
  )
}



# =============================================================================
# 5. plot_species_occupancy_map
# =============================================================================

#' Plot Predicted Occupancy for a Single Species on a Hexgrid Map
#'
#' Produces a quick-check choropleth map of posterior median occupancy
#' probability for one species across the H3 hexgrid. Intended for visual
#' inspection of predictions before saving.
#'
#' @param psi_quantiles Three-dimensional array as returned by
#'   \code{\link{predict_occupancy_chunks}}, with named dimensions
#'   \code{quantile}, \code{species}, and \code{site}.
#' @param species Character. Species name exactly as in
#'   \code{dimnames(psi_quantiles)$species}.
#' @param quantile_name Character. Quantile slice to map.
#'   Defaults to \code{"median"}.
#' @param palette Character. \code{scale_fill_viridis_c} palette option.
#'   Defaults to \code{"rocket"}.
#' @param direction Integer. Viridis palette direction (\code{1} or \code{-1}).
#'   Defaults to \code{-1} (light = high occupancy).
#' @param output_path Character. PNG save path, or \code{NULL} for screen only.
#' @param width,height Numeric. PNG dimensions in cm. Defaults to \code{20x18}.
#' @param dpi Integer. Resolution. Defaults to \code{300}.
#'
#' @return Invisibly returns the \code{ggplot} object.
#'
#' @importFrom h3jsr cell_to_polygon
#' @importFrom ggplot2 ggplot aes geom_sf scale_fill_viridis_c coord_sf labs
#'   theme_minimal theme ggsave element_text element_blank
#' @importFrom ggtext element_markdown
#' @importFrom tibble tibble
#'
#' @export
plot_species_occupancy_map <- function(psi_quantiles,
                                        species,
                                        quantile_name = "median",
                                        palette       = "rocket",
                                        direction     = -1,
                                        output_path   = NULL,
                                        width         = 20,
                                        height        = 18,
                                        dpi           = 300) {

  # ── Validation ────────────────────────────────────────────────────────────────
  stopifnot(
    "species not found in psi_quantiles" =
      species %in% dimnames(psi_quantiles)$species,
    "quantile_name not found in psi_quantiles" =
      quantile_name %in% dimnames(psi_quantiles)$quantile
  )

  # ── Extract values and build hexgrid sf ───────────────────────────────────────
  occ_df <- tibble::tibble(
    h3_address = dimnames(psi_quantiles)$site,
    psi        = psi_quantiles[quantile_name, species, ]
  )

  hex_sf      <- h3jsr::cell_to_polygon(occ_df$h3_address, simple = FALSE)
  hex_sf$psi  <- occ_df$psi

  # ── Build plot ────────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(hex_sf) +
    ggplot2::geom_sf(ggplot2::aes(fill = psi), color = NA) +
    ggplot2::scale_fill_viridis_c(
      option    = palette,
      direction = direction,
      limits    = c(0, 1),
      name      = "Occupancy\nprobability"
    ) +
    ggplot2::coord_sf() +
    ggplot2::labs(
      title    = paste0("Predicted Occupancy \u2014 *", species, "*"),
      subtitle = paste0("Posterior ", quantile_name)
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title      = ggtext::element_markdown(face = "bold"),
      plot.subtitle   = ggplot2::element_text(colour = "grey40", size = 11),
      panel.grid      = ggplot2::element_blank(),
      legend.position = "right"
    )

  # ── Output ────────────────────────────────────────────────────────────────────
  if (!is.null(output_path)) {
    ggplot2::ggsave(output_path, p,
                   width = width, height = height,
                   units = "cm", dpi = dpi, bg = "white")
    message("Map saved: ", basename(output_path))
  } else {
    print(p)
  }

  invisible(p)
}
