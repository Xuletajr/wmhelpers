#' Aggregate Detection Array into Fixed-Width Sampling Windows
#'
#' Collapses the temporal dimension of a 3D presence/absence detection array
#' into fixed-width windows (e.g., 5-day or 7-day periods). Within each
#' window, a species is considered detected (1) if it was detected on at
#' least one day, not detected (0) if sampled but never detected, and
#' missing (NA) if all days in the window were unsampled (NA).
#'
#' The last window may contain fewer days than \code{window} if the total
#' number of days is not an exact multiple of \code{window}. This partial
#' window is retained to avoid discarding valid sampling effort.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @param y A 3D array of dimensions \code{[species × sites × days]}
#'   containing binary detection/non-detection data (0/1) and structural
#'   NAs for unsampled site-days. Typically the output of
#'   \code{\link{build_detection_array}} followed by
#'   \code{\link{build_reindexed_visits}}.
#' @param window Integer. Number of days to aggregate per sampling occasion.
#'   Defaults to \code{5}. If the total number of days is not an exact
#'   multiple of \code{window}, the last window will contain the remaining
#'   days (fewer than \code{window}).
#'
#' @return A 3D array of dimensions \code{[species × sites × n_windows]}
#'   with aggregated detection histories, where
#'   \code{n_windows = ceiling(n_days / window)}.
#'   The last window may span fewer than \code{window} days for sites
#'   whose total sampled days are not an exact multiple of \code{window}.
#'
#' @details
#' All windows are retained, including a partial last window if the total
#' number of days is not divisible by \code{window}. For example, with
#' \code{window = 5} and 23 sampled days, the function produces 5 windows:
#' days 1–5, 6–10, 11–15, 16–20, and 21–23 (partial, 3 days only).
#'
#' This behaviour differs from using \code{floor(n_days / window)}, which
#' would discard the last 3 days entirely. Using \code{ceiling} ensures no
#' valid sampling effort is lost, at the cost of a slightly unbalanced last
#' window for some sites. This is ecologically defensible and transparent in
#' the methods, and \code{spOccupancy::msPGOcc} handles variable effort via
#' structural NAs.
#'
#' A window is assigned \code{NA} for all species if all days within that
#' window are \code{NA} for a given site — meaning the site was not
#' operational during that period.
#'
#' \strong{Important:} \code{\link{build_reindexed_visits}} must be applied
#' \emph{before} this function. Aggregating without prior reindexing inflates
#' the number of windows because leading and trailing NAs from the global
#' time axis are included in the window calculation. See
#' \code{\link{build_reindexed_visits}} for details.
#'
#' Sampling effort (number of valid days per window per site) should be
#' computed after calling this function, using the aggregated array directly:
#'
#' \preformatted{
#' y_reindexed <- build_reindexed_visits(y)
#' y_agg       <- aggregate_detection_windows(y_reindexed, window = 5)
#'
#' # Effort: sites x windows (uses species 1 as NA mask proxy)
#' effort          <- apply(!is.na(y_agg[1, , ]), c(1, 2), as.integer)
#' visits_per_site <- rowSums(effort)
#' }
#'
#' @examples
#' \dontrun{
#' # Correct order: reindex first, then aggregate
#' y_reindexed <- build_reindexed_visits(y)
#' y_agg       <- aggregate_detection_windows(y_reindexed, window = 5)
#' dim(y_agg)
#'
#' # Compute effort
#' effort          <- apply(!is.na(y_agg[1, , ]), c(1, 2), as.integer)
#' visits_per_site <- rowSums(effort)
#' summary(visits_per_site)
#' }
#'
#' @seealso \code{\link{build_detection_array}},
#'   \code{\link{build_reindexed_visits}}
#'
#' @export
aggregate_detection_windows <- function(y, window = 5) {

  n_sp    <- dim(y)[1]
  n_sites <- dim(y)[2]
  n_days  <- dim(y)[3]

  # Use ceiling to include partial last window — no valid days are discarded
  n_windows <- ceiling(n_days / window)

  if (n_windows == 0) {
    stop("No windows fit within ", n_days, " available days. ",
         "Consider reducing 'window'.")
  }

  y_agg <- array(NA, dim = c(n_sp, n_sites, n_windows))

  for (w in seq_len(n_windows)) {

    idx_start <- (w - 1) * window + 1
    idx_end   <- min(w * window, n_days)  # cap at n_days for partial last window
    idx       <- idx_start:idx_end

    # Presence in window = 1 if detected on any day; NA if all days are NA
    y_agg[, , w] <- apply(
      y[, , idx, drop = FALSE], c(1, 2),
      function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    )
  }

  return(y_agg)
}
