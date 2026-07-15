#' Aggregate Detection Array into Fixed-Width Sampling Windows
#'
#' Collapses the temporal dimension of a 3D presence/absence detection array
#' into fixed-width windows (e.g., 5-day or 7-day periods). Within each
#' window, a species is considered detected (1) if it was detected on at
#' least one day, not detected (0) if sampled but never detected, and
#' missing (NA) if all days in the window were unsampled (NA).
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @param y A 3D array of dimensions \code{[species × sites × days]}
#'   containing binary detection/non-detection data (0/1) and structural
#'   NAs for unsampled site-days. Typically the output of
#'   \code{\link{build_detection_array}} followed by
#'   \code{\link{build_reindexed_visits}}.
#' @param window Integer. Number of days to aggregate per sampling occasion.
#'   Defaults to \code{5}. Only complete windows are retained; trailing
#'   days that do not fill a complete window are discarded.
#'
#' @return A 3D array of dimensions \code{[species × sites × n_windows]}
#'   with aggregated detection histories, where \code{n_windows =
#'   floor(n_days / window)}.
#'
#' @details
#' Only complete windows of size \code{window} are retained.
#' For example, with \code{window = 5} and 17 sampled days, the function
#' produces 3 windows (days 1–5, 6–10, 11–15) and discards days 16–17.
#' This behavior is intentional to maintain equal window sizes across sites.
#'
#' A window is assigned \code{NA} for all species if all days within that
#' window are \code{NA} for a given site — meaning the site was not
#' operational during that period.
#'
#' Sampling effort (number of valid days per window per site) should be
#' computed after calling \code{\link{build_reindexed_visits}} on the
#' output of this function, using the reindexed array directly:
#'
#' \preformatted{
#' y_agg <- aggregate_detection_windows(y, window = 5)
#' y     <- build_reindexed_visits(y_agg)
#'
#' # Effort: sites x visits (uses species 1 as NA mask proxy)
#' effort <- apply(!is.na(y[1, , ]), c(1, 2), as.integer)
#' }
#'
#' @examples
#' \dontrun{
#' # Aggregate into 5-day windows
#' y_agg <- aggregate_detection_windows(y_reindexed, window = 5)
#' dim(y_agg)
#'
#' # Then reindex to remove empty windows per site
#' y <- build_reindexed_visits(y_agg)
#'
#' # Compute effort after reindexing
#' effort <- apply(!is.na(y[1, , ]), c(1, 2), as.integer)
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

  # Number of complete windows only — trailing days discarded
  n_windows <- floor(n_days / window)

  if (n_windows == 0) {
    stop("No complete windows of size ", window,
         " fit within the ", n_days, " available days. ",
         "Consider reducing 'window'.")
  }

  y_agg <- array(NA, dim = c(n_sp, n_sites, n_windows))

  for (w in seq_len(n_windows)) {

    idx <- ((w - 1) * window + 1):(w * window)

    # Presence in window = 1 if detected on any day; NA if all days are NA
    y_agg[, , w] <- apply(
      y[, , idx, drop = FALSE], c(1, 2),
      function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    )
  }

  return(y_agg)
}
