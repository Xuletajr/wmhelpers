#' Validate detection array consistency across multiple data sources
#'
#' Validates whether a 3D detection array (`species x sites x dates`)
#' matches reconstructed detections derived from one or more raw
#' monitoring data sources (e.g., acoustic, camera trap, eDNA,
#' citizen science), accounting for deployment periods.
#'
#' This function is intended for quality control prior to occupancy
#' modeling or downstream ecological analyses.
#'
#' @param species_name Character string. Scientific name of the species.
#' @param y_all 3D array. Detection array with dimensions
#'   `species x sites x replicate (date)`.
#' @param date_seq Vector of dates corresponding to the third dimension of `y_all`.
#' @param detection_sources Named list of data frames containing raw detections.
#'   Each data frame must include columns:
#'   `species`, `point_name`, and `date`.
#' @param deployment Data frame containing deployment periods.
#'   Must include: `point_name`, `start_date`, `end_date`.
#' @param verbose Logical. If `TRUE`, prints validation messages per site.
#'
#' @return A tibble summarizing, for each site:
#'   \describe{
#'     \item{species}{Species name}
#'     \item{site}{Site identifier}
#'     \item{mismatches}{Number of mismatched dates}
#'     \item{total_dates}{Total number of evaluated dates}
#'   }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Extracts species-site slices from the detection array
#'   \item Reconstructs expected detections from all provided sources
#'   \item Masks dates outside deployment periods as `NA`
#'   \item Compares reconstructed detections against array values
#' }
#'
#' A mismatch occurs when detection values differ or NA structure differs.
#'
#' @author Thiago Bicudo
#'
#' @export
validate_detection_array <- function(species_name,
                                     y_all,
                                     date_seq,
                                     detection_sources,
                                     deployment,
                                     verbose = TRUE) {

  if (!species_name %in% dimnames(y_all)[[1]]) {
    stop("Species not found in detection array.")
  }

  if (!is.list(detection_sources) || length(detection_sources) == 0) {
    stop("detection_sources must be a non-empty named list of data frames.")
  }

  sites <- dimnames(y_all)[[2]]

  results <- purrr::map_dfr(sites, function(site) {

    y_slice <- tibble::tibble(
      date = date_seq,
      detection_array = as.integer(y_all[species_name, site, ])
    )

    # Reconstruct detections from all sources
    raw_sources <- purrr::map(
      detection_sources,
      ~ .x %>%
        dplyr::filter(species == species_name,
                      point_name == site) %>%
        dplyr::transmute(date, detection = 1L)
    )

    raw <- dplyr::bind_rows(raw_sources) %>%
      dplyr::distinct(date, .keep_all = TRUE) %>%
      dplyr::right_join(y_slice, by = "date") %>%
      tidyr::replace_na(list(detection = 0L)) %>%
      dplyr::arrange(date)

    deploy <- deployment %>%
      dplyr::filter(point_name == site)

    if (nrow(deploy) == 1) {
      raw <- raw %>%
        dplyr::mutate(
          detection_expected = dplyr::if_else(
            date < deploy$start_date |
              date > deploy$end_date,
            NA_integer_,
            detection
          )
        )
    } else {
      raw <- raw %>%
        dplyr::mutate(detection_expected = detection)
    }

    raw <- raw %>%
      dplyr::mutate(diff = detection_array - detection_expected)

    mismatches <- raw %>%
      dplyr::filter(diff != 0 |
                      (is.na(detection_array) != is.na(detection_expected)))

    if (verbose) {
      if (nrow(mismatches) == 0) {
        message("OK | Species: ", species_name,
                " | Site: ", site,
                " | All detections match.")
      } else {
        message("Mismatch | Species: ", species_name,
                " | Site: ", site,
                " | Found ", nrow(mismatches), " mismatches.")
      }
    }

    tibble::tibble(
      species = species_name,
      site = site,
      mismatches = nrow(mismatches),
      total_dates = nrow(raw)
    )
  })

  results
}
