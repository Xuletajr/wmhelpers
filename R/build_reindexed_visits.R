#' Build a Reindexed Visits Array by Removing Empty Time Steps per Site
#'
#' Compacts a species x sites x time detection array by removing time
#' steps (dates) with no sampling effort at a given site, realigning the
#' remaining visits by sampling order. This produces a visit-indexed array
#' (\code{[species x sites x max_visits]}) suitable as input to
#' occupancy models that expect a fixed number of replicate visits per
#' site, rather than a full calendar date sequence.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @param y A 3-dimensional detection array \code{[species x sites x
#'   dates]}, as returned by \code{\link{build_detection_array}} or
#'   \code{\link{build_detection_array_from_events}}. A date column with
#'   \code{NA} for every species at a given site is treated as "no
#'   sampling effort" and is dropped for that site.
#'
#' @return A 3-dimensional array \code{[species x sites x max_visits]},
#'   where \code{max_visits} is the largest number of valid (non-NA)
#'   sampling occasions found across all sites. Sites with fewer valid
#'   visits than \code{max_visits} are padded with \code{NA} in the
#'   remaining visit slots. Sites with zero valid visits (no sampling
#'   effort recorded at all) remain entirely \code{NA} across all visit
#'   slots, rather than causing an error.
#'
#' @details
#' For each site, dates are treated as "no effort" (dropped) only when
#' every species has \code{NA} on that date — this mirrors the array
#' convention where \code{NA} = no sampling effort, \code{0} = effort
#' with no detection, \code{1} = detected.
#'
#' The function is robust to two edge cases that can otherwise cause
#' silent errors: a single-species array (\code{dim(y)[1] == 1}) and
#' sites with zero valid visits. Both are explicitly guarded against by
#' coercing slices back to matrices with known dimensions and by skipping
#' the assignment step for all-NA sites.
#'
#' @seealso \code{\link{build_detection_array}},
#'   \code{\link{build_detection_array_from_events}},
#'   \code{\link{build_reindexed_dates}}
#'
#' @export
build_reindexed_visits <- function(y) {

  stopifnot(
    "y must be a 3D array" = is.array(y) && length(dim(y)) == 3
  )

  n_species <- dim(y)[1]
  n_sites   <- dim(y)[2]
  n_dates   <- dim(y)[3]

  y_new_list <- vector("list", n_sites)

  for (i in seq_len(n_sites)) {

    # Explicitly coerce back to a matrix [species x dates] — y[, i, ] would
    # silently drop to a vector when n_species == 1 or when only one date
    # remains, breaking colSums() below.
    site_matrix <- matrix(y[, i, ], nrow = n_species, ncol = n_dates)

    # A date has effort at this site if at least one species is non-NA
    keep <- which(colSums(!is.na(site_matrix)) > 0)

    y_new_list[[i]] <- matrix(
      site_matrix[, keep, drop = FALSE],
      nrow = n_species,
      ncol = length(keep)
    )
  }

  max_visits <- max(vapply(y_new_list, ncol, integer(1)))

  y_new <- array(NA, dim = c(n_species, n_sites, max_visits))

  for (i in seq_len(n_sites)) {

    n_visits_i <- ncol(y_new_list[[i]])

    # Sites with zero valid visits (no sampling effort recorded at all)
    # are skipped here and simply remain all-NA in y_new — assigning a
    # zero-column slice would otherwise error.
    if (n_visits_i > 0) {
      y_new[, i, seq_len(n_visits_i)] <- y_new_list[[i]]
    }
  }

  y_new
}
