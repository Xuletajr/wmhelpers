#' Set low-effort sampling days to NA in a detection array
#'
#' @description
#' Replaces all species detections with `NA` for site–day combinations identified
#' as having insufficient sampling effort.
#'
#' This function is useful for cleaning acoustic monitoring detection arrays
#' when certain days have too few recordings (e.g., due to recorder failure,
#' late deployment, or early retrieval). All species entries for the affected
#' site-day are replaced with `NA`, preventing those days from being interpreted
#' as true non-detections.
#'
#' @param y_array A 3D array containing detection/nondetection data with
#'   dimensions **species × sites × date**.
#' @param problem_days A data frame listing site–day combinations that should be
#'   treated as missing due to low sampling effort.
#' @param site_col Character string indicating the column in `problem_days`
#'   containing site identifiers matching `dimnames(y_array)$sites`.
#'   Default = `"point_name"`.
#' @param date_col Character string indicating the column in `problem_days`
#'   containing dates matching `dimnames(y_array)$date`.
#'   Default = `"date_array"`.
#'
#' @details
#' For each site–day listed in `problem_days`, the function replaces the entire
#' vector of species detections with `NA`. This ensures that days with
#' insufficient sampling effort are excluded from downstream analyses such as
#' occupancy modeling or detection summaries.
#'
#' @return
#' A modified version of `y_array` where detections for low-effort site–days
#' have been replaced with `NA`.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @export
set_low_effort_days_to_na <- function(y_array,
                                      problem_days,
                                      site_col = "point_name",
                                      date_col = "date_array") {

  # Basic checks
  if(length(dim(y_array)) != 3) {
    stop("y_array must be a 3D array (species × sites × date).")
  }

  if(is.null(dimnames(y_array))) {
    stop("y_array must have dimnames defined.")
  }

  if(!all(c(site_col, date_col) %in% names(problem_days))) {
    stop("Specified site_col or date_col not found in problem_days.")
  }

  # Match indices
  idx_site <- match(problem_days[[site_col]],
                    dimnames(y_array)$sites)

  idx_date <- match(problem_days[[date_col]],
                    dimnames(y_array)$date)

  # Keep valid matches
  valid <- !is.na(idx_site) & !is.na(idx_date)

  # Apply NA replacement
  for(i in which(valid)) {
    y_array[ , idx_site[i], idx_date[i]] <- NA
  }

  return(y_array)
}
