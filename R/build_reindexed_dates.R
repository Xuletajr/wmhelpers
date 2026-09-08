#' Build a Reindexed Sampling Dates Matrix Aligned to Valid Visits
#'
#' Extracts and reindexes sampling dates from a detection array into a
#' \code{[sites x visits]} matrix, dropping time steps with no sampling
#' effort at a given site. The resulting matrix is aligned to the visit
#' structure produced by \code{\link{build_reindexed_visits}} — visit
#' \code{k} in this matrix corresponds to the same sampling occasion as
#' visit \code{k} in the reindexed detection array.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @param y A 3-dimensional detection array \code{[species x sites x
#'   dates]}, as returned by \code{\link{build_detection_array}} or
#'   \code{\link{build_detection_array_from_events}}. Must have non-NULL
#'   \code{dimnames} for the sites (2nd) and dates (3rd) dimensions. A
#'   date column with \code{NA} for every species at a given site is
#'   treated as "no sampling effort" and is dropped for that site.
#'
#' @return A character matrix \code{[sites x max_visits]} of sampling
#'   dates, where \code{max_visits} is the largest number of valid
#'   (non-NA) sampling occasions found across all sites. Sites with fewer
#'   valid visits than \code{max_visits} are padded with \code{NA} in the
#'   remaining visit columns. Sites with zero valid visits (no sampling
#'   effort recorded at all) remain entirely \code{NA} across all visit
#'   columns, rather than causing an error. Row names are site names;
#'   column names are \code{"rep_1"}, \code{"rep_2"}, etc.
#'
#' @details
#' For each site, dates are treated as "no effort" (dropped) only when
#' every species has \code{NA} on that date — this mirrors the array
#' convention where \code{NA} = no sampling effort, \code{0} = effort
#' with no detection, \code{1} = detected. This is the same effort-per-date
#' logic used by \code{\link{build_reindexed_visits}}, so the two
#' functions produce aligned outputs when applied to the same \code{y}.
#'
#' The function is robust to two edge cases that can otherwise cause
#' silent errors: a single-species array (\code{dim(y)[1] == 1}) and
#' sites with zero valid visits. Both are explicitly guarded against by
#' coercing slices back to matrices with known dimensions and by skipping
#' the assignment step for all-NA sites.
#'
#' @seealso \code{\link{build_reindexed_visits}},
#'   \code{\link{build_detection_array}},
#'   \code{\link{build_detection_array_from_events}}
#'
#' @export
build_reindexed_dates <- function(y) {

  stopifnot(
    "y must be a 3D array"                       = is.array(y) &&
      length(dim(y)) == 3,
    "y must have dimnames for the sites dimension" = !is.null(dimnames(y)[[2]]),
    "y must have dimnames for the dates dimension" = !is.null(dimnames(y)[[3]])
  )

  n_species  <- dim(y)[1]
  n_sites    <- dim(y)[2]
  n_dates    <- dim(y)[3]
  site_names <- dimnames(y)[[2]]
  date_names <- dimnames(y)[[3]]

  date_list <- vector("list", n_sites)

  for (i in seq_len(n_sites)) {

    # Explicitly coerce back to a matrix [species x dates] — y[, i, ] would
    # silently drop to a vector when n_species == 1 or when only one date
    # remains, breaking colSums() below.
    site_matrix <- matrix(y[, i, ], nrow = n_species, ncol = n_dates)

    # A date has effort at this site if at least one species is non-NA
    keep <- which(colSums(!is.na(site_matrix)) > 0)

    date_list[[i]] <- date_names[keep]
  }

  max_visits <- max(vapply(date_list, length, integer(1)))

  date_mat <- matrix(NA_character_, nrow = n_sites, ncol = max_visits)

  for (i in seq_len(n_sites)) {

    n_dates_i <- length(date_list[[i]])

    # Sites with zero valid visits (no sampling effort recorded at all)
    # are skipped here and simply remain all-NA in date_mat — assigning a
    # zero-length slice would otherwise error.
    if (n_dates_i > 0) {
      date_mat[i, seq_len(n_dates_i)] <- date_list[[i]]
    }
  }

  rownames(date_mat) <- site_names
  colnames(date_mat) <- paste0("rep_", seq_len(max_visits))

  date_mat
}
