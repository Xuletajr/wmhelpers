#' Read and clean camera trap detections
#'
#' Reads all camera trap CSV files starting with "images"
#' and standardizes species detections.
#'
#' @param config List with project paths.
#' @param region Subfolder inside camtrap directory (e.g., "Alalapadu").
#' @param issue_cam_sites Optional vector of sites to exclude.
#' @param valid_sites Optional vector of valid sites to keep.
#'
#' @return Tibble with columns: sites, species, date
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @export
read_camtrap_detections <- function(config,
                                    region,
                                    issue_cam_sites = NULL,
                                    valid_sites = NULL) {

  cam_dir <- file.path(config$camtrap, region)

  files <- list.files(
    cam_dir,
    pattern = "^images.*\\.csv$",
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop("No camera trap files matching 'images*.csv' found.")
  }

  purrr::map_dfr(files, ~ readr::read_csv(.x, show_col_types = FALSE)) %>%
    dplyr::rename_with(stringr::str_to_lower) %>%
    dplyr::mutate(
      sites = tolower(deployment_id),
      species = paste(genus, species),
      date = as.Date(timestamp)
    ) %>%
    tidyr::drop_na(species) %>%
    dplyr::filter(
      !species %in% c("IDD", "None", "Unknown",
                      "no CV results", "No CV Result")
    ) %>%
    {
      if (!is.null(valid_sites)) {
        dplyr::filter(., sites %in% valid_sites)
      } else .
    } %>%
    {
      if (!is.null(issue_cam_sites)) {
        dplyr::filter(., !sites %in% issue_cam_sites)
      } else .
    } %>%
    dplyr::distinct(sites, date, species, .keep_all = FALSE) %>%
    dplyr::select(sites, species, date) %>%
    dplyr::arrange(sites)
}
