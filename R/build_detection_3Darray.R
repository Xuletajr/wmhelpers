#' Build 3D Detection Array for Multi-Species Occupancy Models
#'
#' Creates a 3D detection array (species × sites × dates) suitable for
#' multi-species occupancy models such as `msPGOcc()` and `intMsPGOcc()`.
#'
#' The array is initialized with `NA` (no sampling effort). Dates within each
#' deployment period are set to `0` (sampling occurred, no detection). Observed
#' detections are set to `1`.
#'
#' Dates containing only `NA` values (i.e., no sampling effort across all sites)
#' are automatically removed from the final array.
#'
#' @param detections A data.frame or tibble containing detection records.
#' @param deployment A data.frame or tibble containing deployment information.
#' Must contain columns: `sites`, `point_name`, `start_date`, `end_date`.
#' @param species_col Character. Column name in `detections` with species IDs.
#' @param site_col Character. Column name in `detections` with site IDs.
#' @param date_col Character. Column name in `detections` with detection dates.
#' @param date_seq A vector of Dates defining replicate sampling occasions.
#'
#' @return A named list with:
#' \describe{
#'   \item{y}{3D integer array (species × sites × dates)}
#'   \item{species}{Character vector of species names}
#'   \item{deployment}{Deployment data.frame used}
#'   \item{date_seq}{Filtered date sequence used in the array}
#' }
#'
#' @details
#' Coding scheme:
#' \itemize{
#'   \item 1 = detection
#'   \item 0 = sampling effort, no detection
#'   \item NA = no sampling effort
#' }
#'
#' @examples
#' \dontrun{
#' audio_array <- build_detection_array(
#'   detections = audioData,
#'   deployment = deploy_audio_filt,
#'   species_col = "species",
#'   site_col = "sites",
#'   date_col = "date",
#'   date_seq = date_seq_aud
#' )
#' }
#'
#' @export
build_detection_array <- function(
    detections,
    deployment,
    species_col,
    site_col,
    date_col,
    date_seq
) {

  stopifnot(
    all(c(species_col, site_col, date_col) %in% names(detections)),
    all(c("sites", "point_name", "start_date", "end_date") %in% names(deployment))
  )

  detections[[date_col]] <- as.Date(detections[[date_col]])
  deployment$start_date <- as.Date(deployment$start_date)
  deployment$end_date   <- as.Date(deployment$end_date)
  date_seq <- as.Date(date_seq)

  species_vec <- sort(unique(detections[[species_col]]))
  N <- length(species_vec)
  J <- nrow(deployment)
  K <- length(date_seq)

  y <- array(
    NA_integer_,
    dim = c(N, J, K),
    dimnames = list(
      species_vec,
      deployment$point_name,
      as.character(date_seq)
    )
  )

  # Mark sampling effort (0) within deployment periods
  for (j in seq_len(J)) {
    active_dates <- date_seq[
      date_seq >= deployment$start_date[j] &
        date_seq <= deployment$end_date[j]
    ]
    date_index <- match(active_dates, date_seq)
    y[, j, date_index] <- 0L
  }

  # Merge detections with deployment
  detections <- merge(
    detections,
    deployment,
    by.x = site_col,
    by.y = "sites"
  )

  detections <- detections[
    detections[[date_col]] >= detections$start_date &
      detections[[date_col]] <= detections$end_date,
  ]

  spp_index  <- match(detections[[species_col]], species_vec)
  site_index <- match(detections$point_name, deployment$point_name)
  date_index <- match(detections[[date_col]], date_seq)

  valid <- !is.na(spp_index) & !is.na(site_index) & !is.na(date_index)

  y[cbind(
    spp_index[valid],
    site_index[valid],
    date_index[valid]
  )] <- 1L

  # Remove dates with only NA
  na_days <- apply(
    y,
    MARGIN = 3,
    FUN = function(x) all(is.na(x))
  )

  if (any(na_days)) {
    y <- y[, , !na_days, drop = FALSE]
    date_seq <- date_seq[!na_days]
  }

  list(
    y = y,
    species = species_vec,
    deployment = deployment,
    date_seq = date_seq
  )
}
