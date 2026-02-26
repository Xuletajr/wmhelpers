#' Read and clean audio detections
#'
#' Reads confirmed audio detections, standardizes site, species and date,
#' filters detections outside the deployment period, and optionally exports
#' detections that occurred outside the valid deployment window.
#'
#' @param config List with project paths. Must contain `species_id`,
#'   `project_name`, and optionally `output_tables`.
#' @param deployment_audio Deployment table with columns:
#'   `sites`, `start_date`, and `end_date`.
#' @param issue_audio_sites Optional character vector of sites to exclude.
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
#' If `export_out_of_period = TRUE`, a CSV file named
#' `"audio_detections_outside_deployment_<project_name>.csv"`
#' will be saved in `output_path`, containing all detections
#' that occurred before `start_date` or after `end_date`.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @export
#' @export
read_audio_detections <- function(config,
                                  deployment_audio,
                                  issue_audio_sites = NULL,
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

  raw_data <- readr::read_csv(file_path,
                              show_col_types = FALSE) %>%
    dplyr::rename(
      species = class_name,
      sites   = site_name
    ) %>%
    dplyr::mutate(
      sites = tolower(sites),
      date  = as.Date(paste(year, month, day, sep = "-"))
    )

  if (!is.null(issue_audio_sites)) {
    raw_data <- dplyr::filter(raw_data,
                              !sites %in% issue_audio_sites)
  }

  raw_data <- raw_data %>%
    dplyr::distinct(sites, date, species)

  joined <- raw_data %>%
    dplyr::left_join(
      deployment_audio %>%
        dplyr::select(sites, start_date, end_date),
      by = "sites"
    )

  # 🔎 Identify out of period
  out_of_period <- joined %>%
    dplyr::filter(date < start_date | date > end_date)

  n_out <- nrow(out_of_period)

  if (export_out_of_period) {

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

    message("📁 Out-of-period detections file exported.")
  }

  # ✅ Keep only valid detections
  joined %>%
    dplyr::filter(date >= start_date &
                    date <= end_date) %>%
    dplyr::select(sites, species, date) %>%
    dplyr::arrange(sites)
}
