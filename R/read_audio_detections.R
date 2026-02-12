#' Read and clean audio detections
#'
#' Reads confirmed audio detections and standardizes site, species and date.
#'
#' @param config List with project paths.
#' @param deployment_audio Deployment table with start and end dates.
#' @param issue_audio_sites Optional vector of sites to exclude.
#'
#' @return Tibble with columns: sites, species, date
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @export
read_audio_detections <- function(config,
                                  deployment_audio,
                                  issue_audio_sites = NULL) {

  file_path <- file.path(
    config$species_id,
    "detections",
    "detections_confirmed.csv"
  )

  readr::read_csv(file_path, show_col_types = FALSE) %>%
    dplyr::rename(
      species = class_name,
      sites   = site_name
    ) %>%
    dplyr::mutate(
      sites = tolower(sites),
      date  = as.Date(paste(year, month, day, sep = "-"))
    ) %>%
    {
      if (!is.null(issue_audio_sites)) {
        dplyr::filter(., !sites %in% issue_audio_sites)
      } else .
    } %>%
    dplyr::distinct(sites, date, species, .keep_all = FALSE) %>%
    dplyr::left_join(
      deployment_audio %>%
        dplyr::select(sites, start_date, end_date),
      by = "sites"
    ) %>%
    dplyr::filter(date >= start_date & date <= end_date) %>%
    dplyr::select(sites, species, date) %>%
    dplyr::arrange(sites)
}
