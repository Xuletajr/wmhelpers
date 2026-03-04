#' Build a 3D Detection Array for Multi-species Occupancy Models
#'
#' Constructs a species × sites × dates detection array from
#' detection records and deployment metadata.
#'
#' The function creates a standardized 3D array where:
#' \itemize{
#'   \item \strong{NA} = no sampling effort (outside deployment period)
#'   \item \strong{0}  = sampling effort occurred but no detection
#'   \item \strong{1}  = species detected
#' }
#'
#' All species share the same sampling-effort structure (NA mask),
#' ensuring consistency required for hierarchical multi-species
#' occupancy models.
#'
#' @param df A data frame of detections containing:
#'   \describe{
#'     \item{species}{Character. Species name.}
#'     \item{sites}{Character. Site identifier.}
#'     \item{date}{Date. Sampling date.}
#'   }
#'
#' @param deploy A data frame containing deployment metadata with:
#'   \describe{
#'     \item{sites}{Character. Site identifier.}
#'     \item{start_date}{Date. Start of sampling period.}
#'     \item{end_date}{Date. End of sampling period.}
#'   }
#'
#' @param site_col Character. Name of the column representing site IDs.
#'   Default is "sites". Can be set to "point_name" if needed.
#'
#' @details
#' The temporal dimension of the array is defined globally using the
#' minimum start date and maximum end date across all deployments.
#'
#' For each site:
#'   \itemize{
#'     \item Dates within the deployment interval are initialized as 0
#'           (sampling effort present).
#'     \item Dates outside the deployment interval remain NA.
#'   }
#'
#' Species detections are then inserted as 1 in the corresponding
#' species-site-date cells.
#'
#' This ensures:
#' \itemize{
#'   \item Identical NA structure across species
#'   \item Site-specific sampling durations
#'   \item Correct 0/1 detection coding
#' }
#'
#' @return A list with:
#'   \describe{
#'     \item{y}{3D array (species × sites × dates).}
#'     \item{species}{Character vector of species names.}
#'     \item{sites}{Character vector of site names.}
#'     \item{date_seq}{Date vector of full temporal sequence.}
#'   }
#'
#' @examples
#' \dontrun{
#' audio_array <- build_detection_array(
#'   df     = audio_filt,
#'   deploy = deploy_audio_filt
#' )
#' }
#'
#' @export

build_detection_array <- function(df,
                                  deploy,
                                  site_col = "sites") {

  # Check if column exists
  if (!site_col %in% names(df)) {
    stop(paste("Column", site_col, "not found in df"))
  }
  if (!site_col %in% names(deploy)) {
    stop(paste("Column", site_col, "not found in deploy"))
  }

  # Standardize types
  df$species <- as.character(df$species)
  df[[site_col]] <- as.character(df[[site_col]])
  df$date <- as.Date(df$date)

  deploy[[site_col]] <- as.character(deploy[[site_col]])
  deploy$start_date <- as.Date(deploy$start_date)
  deploy$end_date <- as.Date(deploy$end_date)

  # Dimensions
  species <- sort(unique(df$species))
  sites   <- sort(unique(deploy[[site_col]]))

  # Global date sequence
  full_dates <- seq(min(deploy$start_date),
                    max(deploy$end_date),
                    by = "day")

  # Initialize array
  y <- array(NA,
             dim = c(length(species),
                     length(sites),
                     length(full_dates)),
             dimnames = list(species = species,
                             sites   = sites,
                             date    = as.character(full_dates)))

  # -----------------
  # 1️⃣ Mark effort
  # -----------------

  idx_site  <- match(deploy[[site_col]], sites)
  idx_start <- match(deploy$start_date, full_dates)
  idx_end   <- match(deploy$end_date, full_dates)

  for (i in seq_len(nrow(deploy))) {
    y[, idx_site[i], idx_start[i]:idx_end[i]] <- 0
  }

  # -----------------
  # 2️⃣ Insert detections
  # -----------------

  idx_sp   <- match(df$species, species)
  idx_site <- match(df[[site_col]], sites)
  idx_date <- match(df$date, full_dates)

  y[cbind(idx_sp, idx_site, idx_date)] <- 1

  return(list(
    y        = y,
    species  = species,
    sites    = sites,
    date_seq = full_dates
  ))
}
