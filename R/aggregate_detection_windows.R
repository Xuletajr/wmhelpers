#' Aggregate Detection Array into Fixed-Width Sampling Windows
#'
#' Collapses the temporal dimension of a 3D presence/absence detection array
#' into fixed-width windows (e.g., 5-day or 7-day periods). Within each
#' window, a species is considered detected (1) if it was detected on at
#' least one day, not detected (0) if sampled but never detected, and
#' missing (NA) if all days in the window were unsampled (NA).
#'
#' This aggregation reduces the size of the detection array, decreases
#' sparsity in detection histories, and typically improves MCMC convergence
#' in multi-species occupancy models (e.g., \code{spOccupancy::msPGOcc}).
#'
#' @param y A 3D array of dimensions \code{[species × sites × days]}
#'   containing binary detection/non-detection data (0/1) and structural
#'   NAs for unsampled site-days. Typically the output of
#'   \code{\link{build_detection_array}} followed by
#'   \code{\link{build_reindexed_visits}}.
#' @param window Integer. Number of days to aggregate per sampling occasion.
#'   Defaults to \code{5}. Only complete windows are retained; trailing
#'   days that do not fill a complete window are discarded.
#' @param export_effort Logical. If \code{TRUE}, returns a list containing
#'   both the aggregated array and a \code{sites × windows} matrix with the
#'   number of valid (non-NA) days per window per site. Useful for
#'   identifying windows with reduced sampling effort, particularly the last
#'   window. Defaults to \code{FALSE}.
#'
#' @return If \code{export_effort = FALSE} (default): a 3D array of
#'   dimensions \code{[species × sites × n_windows]} with aggregated
#'   detection histories.
#'
#'   If \code{export_effort = TRUE}: a named list with two elements:
#'   \describe{
#'     \item{\code{y}}{The aggregated 3D detection array
#'       (\code{species × sites × n_windows}).}
#'     \item{\code{effort}}{A numeric matrix (\code{sites × n_windows})
#'       with the number of valid sampled days per window per site.
#'       Values range from 0 (fully unsampled window) to \code{window}
#'       (fully sampled window). The last window often has fewer days
#'       for sites with non-multiple sampling periods.}
#'   }
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
#' @examples
#' \dontrun{
#' # Basic aggregation into 5-day windows
#' y_agg <- aggregate_detection_windows(y_reindexed, window = 5)
#' dim(y_agg)
#'
#' # With effort matrix to inspect last-window sampling
#' result <- aggregate_detection_windows(y_reindexed, window = 5, export_effort = TRUE)
#' y_agg   <- result$y
#' effort  <- result$effort
#'
#' # Sites where last window has < 5 days
#' effort[, ncol(effort)]
#' }
#'
#' @seealso \code{\link{build_detection_array}},
#'   \code{\link{build_reindexed_visits}}
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @export
aggregate_detection_windows <- function(y, window = 5, export_effort = FALSE) {

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

  y_agg  <- array(NA, dim = c(n_sp, n_sites, n_windows))
  effort <- matrix(0L, nrow = n_sites, ncol = n_windows)

  for (w in seq_len(n_windows)) {

    idx <- ((w - 1) * window + 1):(w * window)

    # Aggregate detection: presence if detected on any day in window
    y_agg[, , w] <- apply(
      y[, , idx, drop = FALSE], c(1, 2),
      function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    )

    # Count valid (non-NA) days per site in this window (using species 1 as mask)
    effort[, w] <- apply(y[1, , idx, drop = FALSE], 1,
                         function(x) sum(!is.na(x)))
  }

  if (export_effort) {
    return(list(y = y_agg, effort = effort))
  }

  return(y_agg)
}
