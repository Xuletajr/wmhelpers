#' Build a Detection Dataset from Dashboard-Derived Detections
#'
#' Filters, standardises, and validates species detections against a
#' deployment table to produce a clean dataset ready for 3D array
#' construction. This function is the dashboard-native counterpart of
#' \code{\link{read_audio_detections}}: instead of reading a file from
#' disk, it receives a detection data frame already loaded and joined
#' in the R environment.
#'
#' The expected input is produced by joining \code{occurrence.csv} with
#' \code{event.csv} and filtering by method and season before calling
#' this function:
#'
#' \preformatted{
#' all_detections <- readr::read_csv(paste0(config$dashboard_export,
#'                                          "occurrence.csv")) %>%
#'   dplyr::filter(occurrenceStatus == "present") %>%
#'   dplyr::rename(species = scientificName) %>%
#'   dplyr::select(eventID, species) %>%
#'   dplyr::left_join(
#'     events %>% dplyr::select(eventID, site, methods, eventDate),
#'     by = "eventID"
#'   )
#'
#' # Season 1 — audio
#' audData_s1 <- wmhelpers::build_detection_data(
#'   detections_df   = all_detections %>%
#'     dplyr::filter(methods == "audio",
#'                   eventDate <= season1_end),
#'   deployment      = deployment_audio_s1,
#'   issue_sites     = issue_audio_sites_s1,
#'   valid_species   = spp_list_audio$species,
#'   project_name    = config$project_name,
#'   method_label    = "audio_s1",
#'   export_out_of_period = TRUE,
#'   output_path     = paste0(config$output_tables, "Summary/")
#' )
#' }
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @param detections_df A data frame of detections already filtered by
#'   method and season. Must contain the columns specified in
#'   \code{site_col}, \code{species_col}, and \code{date_col}. Typically
#'   the result of joining \code{occurrence.csv} with \code{event.csv}
#'   and filtering by \code{methods} and \code{eventDate}.
#' @param deployment A deployment table with columns \code{sites},
#'   \code{start_date}, and \code{end_date}, as returned by
#'   \code{\link{build_deployment_table}} or
#'   \code{\link{read_deployment_file}}.
#' @param site_col Character. Name of the site identifier column in
#'   \code{detections_df}. Site names are converted to lowercase.
#'   Defaults to \code{"site"}.
#' @param species_col Character. Name of the species column in
#'   \code{detections_df}. Defaults to \code{"species"}.
#' @param date_col Character. Name of the date column in
#'   \code{detections_df} (\code{Date} class). Defaults to
#'   \code{"eventDate"}.
#' @param issue_sites Optional character vector of site names to exclude.
#'   Applied after lowercase conversion. Defaults to \code{NULL}.
#' @param species_corrections Optional named character vector for species
#'   name corrections. Names are the original values and values are the
#'   corrected names. Passed to \code{dplyr::recode()}.
#'   Defaults to \code{NULL}.
#' @param valid_species Optional character vector of species names to
#'   retain (e.g., \code{spp_list_audio$species}). If \code{NULL},
#'   all species are kept. Defaults to \code{NULL}.
#' @param project_name Character. Used to name the exported CSV file
#'   when \code{export_out_of_period = TRUE}.
#' @param method_label Character. Short label for the method and season,
#'   used in exported file names (e.g., \code{"audio_s1"},
#'   \code{"camera_s2"}). Defaults to \code{"audio"}.
#' @param export_out_of_period Logical. If \code{TRUE}, detections
#'   outside the deployment window are exported as a CSV. Defaults to
#'   \code{FALSE}.
#' @param output_path Character. Directory where the out-of-period CSV
#'   is saved. Required when \code{export_out_of_period = TRUE}.
#'   Defaults to \code{NULL}.
#'
#' @return A tibble with columns \code{sites}, \code{species}, and
#'   \code{date}, containing only detections that fall within the valid
#'   deployment window of each site. Rows are deduplicated
#'   (one row per unique site × species × date combination) and sorted
#'   by \code{sites}.
#'
#' @details
#' Processing steps:
#' \enumerate{
#'   \item Normalise site names to lowercase.
#'   \item Apply optional site exclusions (\code{issue_sites}).
#'   \item Apply optional species name corrections
#'     (\code{species_corrections}).
#'   \item Apply optional species filter (\code{valid_species}).
#'   \item Deduplicate to unique site × species × date combinations.
#'   \item Join with deployment table to obtain \code{start_date} and
#'     \code{end_date} per site.
#'   \item Identify and optionally export out-of-period detections
#'     (including detections at sites not found in the deployment table).
#'   \item Return only detections within the valid deployment window.
#' }
#'
#' @seealso \code{\link{read_audio_detections}},
#'   \code{\link{build_deployment_table}},
#'   \code{\link{read_deployment_file}},
#'   \code{\link{build_detection_array}}
#'
#' @importFrom dplyr filter mutate rename recode left_join select
#' @importFrom dplyr distinct arrange
#' @importFrom readr write_csv
#'
#' @export
build_detection_data <- function(detections_df,
                                 deployment,
                                 site_col             = "site",
                                 species_col          = "species",
                                 date_col             = "eventDate",
                                 issue_sites          = NULL,
                                 species_corrections  = NULL,
                                 valid_species        = NULL,
                                 project_name         = "",
                                 method_label         = "audio",
                                 export_out_of_period = FALSE,
                                 output_path          = NULL) {

  # ── Input validation ─────────────────────────────────────────────────────────
  stopifnot(
    "detections_df must be a data frame"  = is.data.frame(detections_df),
    "deployment must be a data frame"     = is.data.frame(deployment),
    "site_col must be a single string"    = is.character(site_col) &&
      length(site_col) == 1,
    "species_col must be a single string" = is.character(species_col) &&
      length(species_col) == 1,
    "date_col must be a single string"    = is.character(date_col) &&
      length(date_col) == 1
  )

  required_det  <- c(site_col, species_col, date_col)
  missing_det   <- setdiff(required_det, names(detections_df))
  if (length(missing_det) > 0) {
    stop("Missing columns in detections_df: ",
         paste(missing_det, collapse = ", "))
  }

  required_dep <- c("sites", "start_date", "end_date")
  missing_dep  <- setdiff(required_dep, names(deployment))
  if (length(missing_dep) > 0) {
    stop("Missing columns in deployment: ",
         paste(missing_dep, collapse = ", "))
  }

  # ── Standardise column names and site name case ───────────────────────────
  df <- detections_df %>%
    dplyr::rename(
      sites   = !!rlang::sym(site_col),
      species = !!rlang::sym(species_col),
      date    = !!rlang::sym(date_col)
    ) %>%
    dplyr::mutate(sites = tolower(sites))

  # ── Apply site exclusions ─────────────────────────────────────────────────
  if (!is.null(issue_sites)) {
    df <- df %>% dplyr::filter(!sites %in% tolower(issue_sites))
  }

  # ── Apply species name corrections ────────────────────────────────────────
  if (!is.null(species_corrections)) {
    df <- df %>%
      dplyr::mutate(species = dplyr::recode(species, !!!species_corrections))
  }

  # ── Apply species filter ──────────────────────────────────────────────────
  if (!is.null(valid_species)) {
    df <- df %>% dplyr::filter(species %in% valid_species)
  }

  # ── Deduplicate ───────────────────────────────────────────────────────────
  df <- df %>% dplyr::distinct(sites, species, date)

  # ── Join with deployment to validate period ───────────────────────────────
  joined <- df %>%
    dplyr::left_join(
      deployment %>% dplyr::select(sites, start_date, end_date),
      by = "sites"
    )

  # ── Identify out-of-period detections ────────────────────────────────────
  out_of_period <- joined %>%
    dplyr::filter(
      is.na(start_date) |
        date < start_date |
        date > end_date
    )

  n_out <- nrow(out_of_period)

  if (n_out == 0) {
    message("[", method_label, "] No detections outside deployment period.")
  } else {
    message("[", method_label, "] ", n_out,
            " detection(s) outside deployment period identified.")
  }

  # ── Export out-of-period detections ──────────────────────────────────────
  if (export_out_of_period && n_out > 0) {

    if (is.null(output_path)) {
      stop("`output_path` must be provided when `export_out_of_period = TRUE`.")
    }

    out_file <- file.path(
      output_path,
      paste0(method_label, "_detections_outside_deployment_",
             project_name, ".csv")
    )

    readr::write_csv(out_of_period, file = out_file)
    message("[", method_label, "] Out-of-period detections exported to:\n  ",
            out_file)
  }

  # ── Return valid detections only ──────────────────────────────────────────
  joined %>%
    dplyr::filter(
      !is.na(start_date) &
        date >= start_date &
        date <= end_date
    ) %>%
    dplyr::select(sites, species, date) %>%
    dplyr::arrange(sites)
}
