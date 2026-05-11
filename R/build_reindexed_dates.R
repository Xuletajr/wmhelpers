#' Build a reindexed sampling dates matrix aligned to valid visits
#'
#' Extracts and reindexes sampling dates from a species x sites x time
#' array, returning a 2D matrix (sites x visits) where each row contains
#' only the dates of valid sampling occasions for that site. Valid occasions
#' are those where at least one species has a non-\code{NA} detection value.
#' The result is aligned to the same visit structure produced by
#' \code{build_reindexed_visits()}, making it suitable for use as a
#' date covariate in multi-species occupancy models (e.g., via
#' \code{spOccupancy}).
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @param y A 3-dimensional array of dimensions \code{[species x sites x time]}
#'   with named dimensions. \code{dimnames(y)[[2]]} must contain site names and
#'   \code{dimnames(y)[[3]]} must contain date identifiers (e.g., Julian day
#'   as character strings).
#'
#' @return A character matrix of dimensions \code{[sites x max_visits]}, where
#'   \code{max_visits} is the maximum number of valid visits across all sites.
#'   Row names correspond to site names and column names follow the pattern
#'   \code{rep_1}, \code{rep_2}, ..., \code{rep_K}. Sites with fewer valid
#'   visits are padded with \code{NA}.
#'
#' @details
#' This function is the date-alignment counterpart of
#' \code{build_reindexed_visits()} and must be applied to the same input
#' array \code{y} to ensure that detection histories and date covariates
#' share an identical visit structure. Misalignment between the two outputs
#' will cause errors or silent inconsistencies in downstream occupancy models.
#'
#' A sampling occasion is considered valid for a given site if at least one
#' species has a non-\code{NA} value in that time step of the site matrix.
#'
#' @export
build_reindexed_dates <- function(y) {

  stopifnot(
    is.array(y),
    length(dim(y)) == 3,
    !is.null(dimnames(y)[[2]]),
    !is.null(dimnames(y)[[3]])
  )

  site_names <- dimnames(y)[[2]]
  date_names <- dimnames(y)[[3]]

  sites <- dim(y)[2]

  date_list <- list()

  for(i in 1:sites){

    site_matrix <- y[, i, ]

    # Identify valid sampling occasions (at least one non-NA)
    keep <- which(colSums(!is.na(site_matrix)) > 0)

    # Store dates for this site
    date_list[[i]] <- date_names[keep]
  }

  max_visits <- max(sapply(date_list, length))

  # Initialize matrix
  date_mat <- matrix(NA, nrow = sites, ncol = max_visits)

  for(i in 1:sites){
    date_mat[i, 1:length(date_list[[i]])] <- date_list[[i]]
  }

  rownames(date_mat) <- site_names
  colnames(date_mat) <- paste0("rep_", seq_len(max_visits))

  return(date_mat)
}
