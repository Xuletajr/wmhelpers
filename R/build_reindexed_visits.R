#' Build a reindexed visits array by removing empty time steps per site
#'
#' Converts a species × sites × time array with a global time structure
#' (e.g., calendar dates) into a compacted array where each site retains
#' only its valid sampling events. Empty time steps (i.e., columns where
#' all species have \code{NA}) are removed per site, and the resulting
#' visits are re-aligned by sampling order. The output array is padded
#' with \code{NA} to accommodate the site with the highest number of
#' valid visits.
#'
#' @param y A 3-dimensional array of dimensions \code{[species × sites × time]},
#'   where \code{NA} indicates the absence of a sampling event for a given
#'   species, site, and time step.
#'
#' @return A 3-dimensional array of dimensions \code{[species × sites × max_visits]},
#'   where \code{max_visits} is the maximum number of valid visits recorded
#'   across all sites. Sites with fewer visits are padded with \code{NA}.
#'
#' @details
#' This function is particularly useful when detection data are structured
#' around a global calendar (e.g., daily or weekly time steps) but individual
#' sites were sampled on different and irregular schedules. Removing empty
#' time steps reduces array size and improves computational efficiency in
#' downstream occupancy models, while fully preserving each site's detection
#' history.
#'
#' A time step is considered valid for a given site if at least one species
#' has a non-\code{NA} value in that column of the site matrix.
#'
#' @export
build_reindexed_visits <- function(y) {
  sp <- dim(y)[1]
  sites <- dim(y)[2]

  y_new_list <- list()

  for(i in 1:sites){

    site_matrix <- y[, i, ]

    keep <- which(colSums(!is.na(site_matrix)) > 0)

    y_new_list[[i]] <- site_matrix[, keep]
  }

  max_visits <- max(sapply(y_new_list, ncol))

  y_new <- array(NA, dim = c(sp, sites, max_visits))

  for(i in 1:sites){
    y_new[, i, 1:ncol(y_new_list[[i]])] <- y_new_list[[i]]
  }

  return(y_new)
}
