# =============================================================================
# Detection summary, visualization and export functions
# File: R/detection_summary.R
# =============================================================================
#'
#' Summarise Species Detections
#'
#' Groups a detection data frame by species (and optional grouping variables)
#' and returns the total number of detection records and the number of unique
#' sites per species. Joins a species list to add common name and IUCN status
#' code.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @param detections_df A data frame of detections. Must contain at minimum
#'   a species column (matching \code{species_col}) and a sites column named
#'   \code{"sites"}. Additional columns referenced in \code{grouping_vars}
#'   (e.g., \code{method}, \code{season}, \code{community}, \code{landcover})
#'   must be present if used.
#' @param spp_list A data frame containing species metadata. Must contain at
#'   minimum a column matching \code{species_col} and the columns specified in
#'   \code{common_name_col} and \code{iucn_col}.
#' @param species_col Character. Name of the species column in both
#'   \code{detections_df} and \code{spp_list}. Defaults to \code{"species"}.
#' @param common_name_col Character. Name of the common name column in
#'   \code{spp_list}. Defaults to \code{"common_name"}.
#' @param iucn_col Character. Name of the IUCN status column in
#'   \code{spp_list}. Defaults to \code{"iucn_code"}.
#' @param grouping_vars Optional character vector of additional column names
#'   in \code{detections_df} to group by alongside species (e.g.,
#'   \code{"method"}, \code{c("method", "landcover")},
#'   \code{"community"}). If \code{NULL} (default), results are summarised
#'   only by species.
#'
#' @return A tibble with one row per species (x grouping combination if
#'   \code{grouping_vars} is set) and columns:
#' \describe{
#'   \item{\code{species}}{Scientific name.}
#'   \item{\code{common_name}}{Common name from \code{spp_list}.}
#'   \item{\code{iucn_code}}{IUCN status code from \code{spp_list}.}
#'   \item{\code{<grouping_vars>}}{Additional grouping columns, if provided.}
#'   \item{\code{N_detections}}{Total number of detection records.}
#'   \item{\code{N_sites}}{Number of unique sites where the species
#'     was detected.}
#' }
#'
#' @details
#' This function does not apply any filtering internally — all subsetting
#' (by method, season, community, landcover, etc.) should be done on
#' \code{detections_df} before calling this function, or passed via
#' \code{grouping_vars} to split results by those dimensions.
#'
#' Species present in \code{detections_df} but absent from \code{spp_list}
#' will have \code{NA} in \code{common_name} and \code{iucn_code}.
#'
#' The IUCN column is always renamed to \code{"iucn_code"} in the output,
#' regardless of the value passed to \code{iucn_col}. This ensures consistent
#' column naming for downstream functions such as
#' \code{\link{plot_detection_summary}} and \code{\link{export_detection_plots}}.
#'
#' @examples
#' \dontrun{
#' # Overall by species only
#' summarise_detections(detections, spp_list, iucn_col = "iucn")
#'
#' # By species + method
#' summarise_detections(detections, spp_list, iucn_col = "iucn",
#'                      grouping_vars = "method")
#'
#' # By species + method + landcover
#' summarise_detections(detections, spp_list, iucn_col = "iucn",
#'                      grouping_vars = c("method", "landcover"))
#' }
#'
#' @seealso \code{\link{plot_detection_summary}},
#'   \code{\link{export_detection_plots}}
#'
#' @importFrom dplyr left_join select group_by summarise n n_distinct across all_of rename
#'
#' @export
summarise_detections <- function(detections_df,
                                 spp_list,
                                 species_col     = "species",
                                 common_name_col = "common_name",
                                 iucn_col        = "iucn_code",
                                 grouping_vars   = NULL) {

  # ── Input validation ─────────────────────────────────────────────────────────
  stopifnot(
    "detections_df must be a data frame"  = is.data.frame(detections_df),
    "spp_list must be a data frame"       = is.data.frame(spp_list),
    "species_col must be a single string" = is.character(species_col) &&
      length(species_col) == 1,
    "'sites' column must exist in detections_df" =
      "sites" %in% names(detections_df)
  )

  if (!species_col %in% names(detections_df)) {
    stop("Column '", species_col, "' not found in detections_df.")
  }

  missing_spp <- setdiff(c(species_col, common_name_col, iucn_col),
                         names(spp_list))
  if (length(missing_spp) > 0) {
    stop("Missing columns in spp_list: ",
         paste(missing_spp, collapse = ", "))
  }

  if (!is.null(grouping_vars)) {
    missing_grp <- setdiff(grouping_vars, names(detections_df))
    if (length(missing_grp) > 0) {
      stop("grouping_vars columns not found in detections_df: ",
           paste(missing_grp, collapse = ", "))
    }
  }

  # ── Build group_by keys: species + metadata + optional grouping_vars ─────────
  group_keys <- c(species_col, common_name_col, iucn_col, grouping_vars)

  # ── Join metadata and summarise ───────────────────────────────────────────────
  detections_df %>%
    dplyr::left_join(
      spp_list %>%
        dplyr::select(dplyr::all_of(c(species_col, common_name_col, iucn_col))),
      by = species_col
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_keys))) %>%
    dplyr::summarise(
      N_detections = dplyr::n(),
      N_sites      = dplyr::n_distinct(sites),
      .groups      = "drop"
    ) %>%
    dplyr::rename(
      species     = !!species_col,
      common_name = !!common_name_col,
      iucn_code   = !!iucn_col        # always standardized to "iucn_code"
    )
}


#' Plot Species Detection Summary
#'
#' Produces a pair of bar plots from a species detection summary — one for
#' number of detections and one for number of sites — by calling
#' \code{\link{plot_bar_species}} for each. Returns a named list so the user
#' can inspect or customise each plot before exporting.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @param det_summary A tibble as returned by \code{\link{summarise_detections}}.
#'   Must contain columns \code{common_name}, \code{iucn_code},
#'   \code{N_detections}, and \code{N_sites}.
#' @param top_n Integer. Number of top species to display in each plot.
#'   Defaults to \code{10}.
#' @param color_det Character. Bar fill color for the detections plot.
#'   Defaults to \code{"#D99E73"}.
#' @param color_sites Character. Bar fill color for the sites plot.
#'   Defaults to \code{"#008AA1"}.
#' @param subtitle Character or \code{NULL}. Subtitle applied to both plots,
#'   typically the filter label (e.g., \code{"audio_s1"},
#'   \code{"camera_wet"}). Defaults to \code{NULL}.
#' @param status_col Character. Column name for IUCN status in
#'   \code{det_summary}. Defaults to \code{"iucn_code"}.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\code{detections}}{A \code{ggplot2} object showing the top species
#'     by number of detection records.}
#'   \item{\code{sites}}{A \code{ggplot2} object showing the top species
#'     by number of unique sites.}
#' }
#'
#' @seealso \code{\link{summarise_detections}},
#'   \code{\link{export_detection_plots}},
#'   \code{\link{plot_bar_species}}
#'
#' @export
plot_detection_summary <- function(det_summary,
                                   top_n       = 10,
                                   color_det   = "#D99E73",
                                   color_sites = "#008AA1",
                                   subtitle    = NULL,
                                   status_col  = "iucn_code") {

  # ── Input validation ─────────────────────────────────────────────────────────
  stopifnot(
    "det_summary must be a data frame" = is.data.frame(det_summary),
    "top_n must be a positive integer" = is.numeric(top_n) && top_n > 0
  )

  required_cols <- c("common_name", status_col, "N_detections", "N_sites")
  missing_cols  <- setdiff(required_cols, names(det_summary))
  if (length(missing_cols) > 0) {
    stop("Missing columns in det_summary: ",
         paste(missing_cols, collapse = ", "))
  }

  # ── Detections plot ───────────────────────────────────────────────────────────
  p_det <- wmhelpers::plot_bar_species(
    dat       = det_summary,
    value     = "N_detections",
    species   = "common_name",
    status    = status_col,
    top_n     = top_n,
    bar_color = color_det,
    xlab      = "Number of detections",
    subtitle  = subtitle
  )

  # ── Sites plot ────────────────────────────────────────────────────────────────
  p_sites <- wmhelpers::plot_bar_species(
    dat       = det_summary,
    value     = "N_sites",
    species   = "common_name",
    status    = status_col,
    top_n     = top_n,
    bar_color = color_sites,
    xlab      = "Number of sites",
    subtitle  = subtitle
  )

  list(detections = p_det, sites = p_sites)
}


#' Export Detection Summary Plots
#'
#' High-level pipeline function that iterates over a list of filter
#' specifications, applies each filter to the detection data, calls
#' \code{\link{summarise_detections}} and \code{\link{plot_detection_summary}},
#' and exports each plot as a separate PNG file.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @param detections A data frame of detections containing at minimum
#'   \code{sites}, \code{species}, and any columns referenced in
#'   \code{plot_specs} filters (e.g., \code{method}, \code{season},
#'   \code{community}, \code{landcover}).
#' @param spp_list A species metadata data frame passed to
#'   \code{\link{summarise_detections}}.
#' @param plot_specs A list of plot specifications. Each element must be a
#'   list with two elements:
#'   \describe{
#'     \item{\code{filters}}{A named list of column-value pairs used to filter
#'       \code{detections} before summarising. Use \code{NULL} for no filtering.}
#'     \item{\code{label}}{Character. A short descriptor used in the subtitle
#'       and output filename (e.g., \code{"audio_wet"}).}
#'   }
#' @param output_path Character. Directory where PNG files are saved.
#' @param project_name Character. Appended to output filenames.
#' @param dpi Numeric. Output resolution. Defaults to \code{600}.
#' @param color_det Character. Bar fill color for detections plot.
#'   Defaults to \code{"#D99E73"}.
#' @param color_sites Character. Bar fill color for sites plot.
#'   Defaults to \code{"#008AA1"}.
#' @param top_n Integer. Number of top species per plot. Defaults to \code{10}.
#' @param width Numeric. Output image width in cm. Defaults to \code{18}.
#' @param height Numeric. Output image height in cm. Defaults to \code{16}.
#' @param species_col Character. Species column name. Defaults to
#'   \code{"species"}.
#' @param common_name_col Character. Common name column name in
#'   \code{spp_list}. Defaults to \code{"common_name"}.
#' @param iucn_col Character. IUCN status column name in \code{spp_list}.
#'   Defaults to \code{"iucn_code"}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   writing PNG files to \code{output_path}.
#'
#' @details
#' Each spec produces two separate PNG files:
#' \itemize{
#'   \item \code{species_detection_det_<label>_<project_name>.png}
#'   \item \code{species_detection_sites_<label>_<project_name>.png}
#' }
#'
#' If a filter combination results in zero rows, the iteration is skipped
#' with a warning and processing continues with the next specification.
#'
#' @seealso \code{\link{summarise_detections}},
#'   \code{\link{plot_detection_summary}}
#'
#' @importFrom purrr walk
#' @importFrom grDevices png dev.off
#'
#' @export
export_detection_plots <- function(detections,
                                   spp_list,
                                   plot_specs,
                                   output_path,
                                   project_name,
                                   dpi             = 600,
                                   color_det       = "#D99E73",
                                   color_sites     = "#008AA1",
                                   top_n           = 10,
                                   width           = 18,
                                   height          = 16,
                                   species_col     = "species",
                                   common_name_col = "common_name",
                                   iucn_col        = "iucn_code") {

  # ── Input validation ─────────────────────────────────────────────────────────
  stopifnot(
    "detections must be a data frame"      = is.data.frame(detections),
    "spp_list must be a data frame"        = is.data.frame(spp_list),
    "plot_specs must be a list"            = is.list(plot_specs),
    "output_path must be a single string"  = is.character(output_path) &&
      length(output_path) == 1,
    "project_name must be a single string" = is.character(project_name) &&
      length(project_name) == 1,
    "dpi must be a positive number"        = is.numeric(dpi) && dpi > 0
  )

  if (!endsWith(output_path, "/")) output_path <- paste0(output_path, "/")

  # ── Validate plot_specs structure ─────────────────────────────────────────────
  for (i in seq_along(plot_specs)) {
    spec <- plot_specs[[i]]
    if (!all(c("filters", "label") %in% names(spec))) {
      stop("plot_specs[[", i, "]] must have 'filters' and 'label' elements.")
    }
  }

  # ── Iterate and export ────────────────────────────────────────────────────────
  purrr::walk(plot_specs, function(spec) {

    # ── Apply all filters dynamically ──────────────────────────────────────────
    dat <- detections

    if (!is.null(spec$filters) && length(spec$filters) > 0) {
      for (col in names(spec$filters)) {
        val <- spec$filters[[col]]
        if (!col %in% names(dat)) {
          warning("Filter column '", col, "' not found in detections. ",
                  "Skipping spec '", spec$label, "'.")
          return(invisible(NULL))
        }
        dat <- dplyr::filter(dat, .data[[col]] == val)
      }
    }

    if (nrow(dat) == 0) {
      warning("No rows remaining after filtering for spec '",
              spec$label, "'. Skipping.")
      return(invisible(NULL))
    }

    # ── Summarise ──────────────────────────────────────────────────────────────
    det_summary <- summarise_detections(
      detections_df   = dat,
      spp_list        = spp_list,
      species_col     = species_col,
      common_name_col = common_name_col,
      iucn_col        = iucn_col
    )

    # ── Build plots ────────────────────────────────────────────────────────────
    plots <- plot_detection_summary(
      det_summary = det_summary,
      top_n       = top_n,
      color_det   = color_det,
      color_sites = color_sites,
      subtitle    = spec$label,
      status_col  = "iucn_code"  # always standardized by summarise_detections
    )

    # ── Export detections plot ─────────────────────────────────────────────────
    out_det <- paste0(
      output_path,
      "species_detection_det_", spec$label, "_", project_name, ".png"
    )
    grDevices::png(filename = out_det, width = width, height = height,
                   units = "cm", bg = "white", res = dpi)
    print(plots$detections)
    grDevices::dev.off()
    message("Exported: ", basename(out_det))

    # ── Export sites plot ──────────────────────────────────────────────────────
    out_sites <- paste0(
      output_path,
      "species_detection_sites_", spec$label, "_", project_name, ".png"
    )
    grDevices::png(filename = out_sites, width = width, height = height,
                   units = "cm", bg = "white", res = dpi)
    print(plots$sites)
    grDevices::dev.off()
    message("Exported: ", basename(out_sites))
  })

  invisible(NULL)
}
