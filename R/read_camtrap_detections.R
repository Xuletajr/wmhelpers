#' Read and clean camera trap detections
#'
#' Recursively reads all camera trap CSV files matching
#' `"images*.csv"` inside the camtrap directory, standardizes detections,
#' filters detections outside the deployment period, and optionally exports
#' detections occurring outside the valid deployment window.
#'
#' @param config List with project paths. Must contain `camtrap`,
#'   `project_name`, and optionally `output_tables`.
#' @param deployment_camera Deployment table with columns:
#'   `sites`, `start_date`, and `end_date`.
#' @param issue_cam_sites Optional character vector of sites to exclude.
#' @param valid_sites Optional character vector of valid sites to retain.
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
#' The function:
#' \enumerate{
#'   \item Recursively locates all files matching `"images*.csv"`
#'   \item Standardizes column names and species formatting
#'   \item Removes blank and unidentified records
#'   \item Filters issue or invalid sites if provided
#'   \item Filters detections outside deployment dates
#' }
#'
#' If `export_out_of_period = TRUE`, a CSV file named
#' `"camera_detections_outside_deployment_<project_name>.csv"`
#' will be saved in `output_path`, containing all detections
#' that occurred before `start_date` or after `end_date`.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @export
#' @export
read_camtrap_detections <- function(config,
                                    deployment_camera,
                                    issue_cam_sites = NULL,
                                    valid_sites = NULL,
                                    export_out_of_period = FALSE,
                                    output_path = NULL) {

  cam_root <- normalizePath(config$camtrap,
                            mustWork = TRUE)

  cam_files <- fs::dir_ls(
    path = cam_root,
    recurse = TRUE,
    regexp = "images.*\\.csv$",
    type = "file"
  )

  if (length(cam_files) == 0) {
    stop("No camera trap files matching 'images*.csv' found.")
  }

  raw_data <- purrr::map_dfr(
    cam_files,
    readr::read_csv,
    show_col_types = FALSE
  ) %>%
    dplyr::rename_with(stringr::str_to_lower) %>%
    dplyr::mutate(
      sites   = tolower(deployment_id),
      species = paste(genus, species),
      date    = as.Date(timestamp)
    ) %>%
    tidyr::drop_na(species) %>%
    dplyr::filter(
      !species %in% c("IDD", "None", "Unknown",
                      "no CV results", "No CV Result")
    )

  if (!is.null(valid_sites)) {
    raw_data <- dplyr::filter(raw_data,
                              sites %in% valid_sites)
  }

  if (!is.null(issue_cam_sites)) {
    raw_data <- dplyr::filter(raw_data,
                              !sites %in% issue_cam_sites)
  }

  raw_data <- raw_data %>%
    dplyr::distinct(sites, date, species)

  # 🔎 Join with deployment to validate period
  joined <- detections %>%
    dplyr::left_join(
      deployment_camera %>%
        dplyr::select(sites, start_date, end_date),
      by = "sites"
    )

  # 🔎 Identify out of period (including missing deployment)
  out_of_period <- joined %>%
    dplyr::filter(
      is.na(start_date) |
        date < start_date |
        date > end_date
    )

  n_out <- nrow(out_of_period)

  if (n_out == 0) {
    message("[Camera] ✅ No detections outside deployment period found.")
  } else {
    message("[Camera] ⚠ ", n_out,
            " detections outside deployment period identified.")
  }

  if (export_out_of_period) {

    if (is.null(output_path)) {
      stop("`output_path` must be provided when `export_out_of_period = TRUE`.")
    }

    readr::write_csv(
      out_of_period,
      file = file.path(
        output_path,
        paste0("camera_detections_outside_deployment_",
               config$project_name, ".csv")
      )
    )

    message("[Camera] 📁 Out-of-period detections file exported.")
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
