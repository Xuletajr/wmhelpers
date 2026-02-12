#' Recursively reads all camera trap CSV files matching
#' "images*.csv" inside the camtrap directory and
#' standardizes detections.
#'
#' @param config List with project paths.
#' @param issue_cam_sites Optional vector of sites to exclude.
#' @param valid_sites Optional vector of valid sites to retain.
#'
#' @return Tibble with columns: sites, species, date
#' #' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @export
read_camtrap_detections <- function(config,
                                    issue_cam_sites = NULL,
                                    valid_sites = NULL) {

  # 1️⃣ List all images*.csv files recursively
  cam_files <- fs::dir_ls(
    path = config$camtrap,
    recurse = TRUE,
    regexp = "^images.*\\.csv$",
    type = "file"
  )

  if (length(cam_files) == 0) {
    stop("No camera trap files matching 'images*.csv' found.")
  }

  # 2️⃣ Read and bind all files
  purrr::map_dfr(
    cam_files,
    readr::read_csv,
    show_col_types = FALSE
  ) %>%

    # 3️⃣ Standardize names
    dplyr::rename_with(stringr::str_to_lower) %>%

    # 4️⃣ Create standardized columns
    dplyr::mutate(
      sites   = tolower(deployment_id),
      species = paste(genus, species),
      date    = as.Date(timestamp)
    ) %>%

    # 5️⃣ Remove blanks and unidentified species
    tidyr::drop_na(species) %>%
    dplyr::filter(
      !species %in% c("IDD", "None", "Unknown",
                      "no CV results", "No CV Result")
    ) %>%

    # 6️⃣ Keep only valid sites (if provided)
    {
      if (!is.null(valid_sites)) {
        dplyr::filter(., sites %in% valid_sites)
      } else .
    } %>%

    # 7️⃣ Remove issue sites (if provided)
    {
      if (!is.null(issue_cam_sites)) {
        dplyr::filter(., !sites %in% issue_cam_sites)
      } else .
    } %>%

    # 8️⃣ Keep unique species-site-day
    dplyr::distinct(sites, date, species, .keep_all = FALSE) %>%

    # 9️⃣ Final clean output
    dplyr::select(sites, species, date) %>%
    dplyr::arrange(sites)
}
