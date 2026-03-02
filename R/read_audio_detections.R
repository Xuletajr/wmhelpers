#' Read and clean audio detections
#'
#' Reads confirmed audio detections, standardizes site, species and date,
#' applies optional species corrections and species filtering,
#' filters detections outside the deployment period, and optionally exports
#' detections occurring outside the valid deployment window.
#'
#' @param config List with project paths. Must contain `species_id`
#'   and `project_name`.
#' @param deployment_audio Deployment table with columns:
#'   `sites`, `start_date`, and `end_date`.
#' @param issue_audio_sites Optional character vector of sites to exclude.
#' @param species_corrections Optional named character vector for species
#'   name corrections. Names must be the original values and elements the
#'   corrected names. Passed to `dplyr::recode()`.
#' @param valid_species Optional character vector of species to retain
#'   (e.g., `spp_list_audio$species`).
#' @param export_out_of_period Logical. If `TRUE`, detections occurring
#'   outside the deployment period are exported as a CSV file.
#'   Default is `FALSE`.
#' @param output_path Character path where the out-of-period detections
#'   CSV file will be saved. Required if `export_out_of_period = TRUE`.
#'
#' @return A tibble with columns:
#'   \itemize{
#'     \item `sites`
#'     \item `species`
#'     \item `date`
#'   }
#'   containing only detections within the valid deployment period.
#'
#' @details
#' Processing steps:
#' \enumerate{
#'   \item Reads confirmed audio detections
#'   \item Standardizes column names and date format
#'   \item Applies optional species corrections
#'   \item Applies optional site and species filters
#'   \item Filters detections outside deployment dates
#' }
#'
#' If `export_out_of_period = TRUE` and detections outside the deployment
#' window are found, a CSV file named
#' `"audio_detections_outside_deployment_<project_name>.csv"`
#' will be saved in `output_path`.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @export
read_audio_detections <- function(config,
                                  deployment_audio,
                                  issue_audio_sites = NULL,
                                  species_corrections = NULL,
                                  valid_species = NULL,
                                  export_out_of_period = FALSE,
                                  output_path = NULL) {

  file_path <- normalizePath(
    file.path(
      config$species_id,
      "detections",
      "detections_confirmed.csv"
    ),
    mustWork = TRUE
  )

  # 🔹 Read and standardize raw detections
  raw_data <- readr::read_csv(
    file_path,
    show_col_types = FALSE
  ) %>%
    dplyr::rename(
      species = class_name,
      sites   = site_name
    ) %>%
    dplyr::mutate(
      sites = tolower(sites),
      date  = as.Date(paste(year, month, day, sep = "-"))
    )

  # 🔹 Apply species corrections (if provided)
  if (!is.null(species_corrections)) {
    raw_data <- raw_data %>%
      dplyr::mutate(
        species = dplyr::recode(species, !!!species_corrections)
      )
  }

  # 🔹 Optional site filters
  if (!is.null(issue_audio_sites)) {
    raw_data <- dplyr::filter(raw_data, !sites %in% issue_audio_sites)
  }

  # 🔹 Optional species filter
  if (!is.null(valid_species)) {
    raw_data <- dplyr::filter(raw_data, species %in% valid_species)
  }

  raw_data <- raw_data %>%
    dplyr::distinct(sites, date, species)

  # 🔎 Join with deployment to validate period
  joined <- raw_data %>%
    dplyr::left_join(
      deployment_audio %>%
        dplyr::select(sites, start_date, end_date),
      by = "sites"
    )

  # 🔎 Identify out-of-period (including missing deployment)
  out_of_period <- joined %>%
    dplyr::filter(
      is.na(start_date) |
        date < start_date |
        date > end_date
    )

  n_out <- nrow(out_of_period)

  if (n_out == 0) {
    message("[Audio] ✅ No detections outside deployment period found.")
  } else {
    message("[Audio] ⚠ ", n_out,
            " detections outside deployment period identified.")
  }

  # 📁 Export ONLY if requested AND if detections exist
  if (export_out_of_period && n_out > 0) {

    if (is.null(output_path)) {
      stop("`output_path` must be provided when `export_out_of_period = TRUE`.")
    }

    readr::write_csv(
      out_of_period,
      file = file.path(
        output_path,
        paste0("audio_detections_outside_deployment_",
               config$project_name, ".csv")
      )
    )

    message("[Audio] 📁 Out-of-period detections file exported.")
  }

  # ✅ Keep only valid detections
  joined %>%
    dplyr::filter(
      !is.na(start_date) &
        date >= start_date &
        date <= end_date
    ) %>%
    dplyr::select(sites, species, date) %>%
    dplyr::arrange(sites)
}
