#' Build a 3D Detection Array from Raw Sampling Events
#'
#' Constructs a species x sites x dates detection array from detection
#' records and a raw table of sampling events (ground-truth effort).
#' Unlike \code{\link{build_detection_array}}, which derives the sampled
#' date range from \code{deploy$start_date}/\code{deploy$end_date} and
#' therefore assumes continuous deployment per site, this version uses the
#' raw Darwin-Core-style events table as the source of sampling effort —
#' one row per actual recording or sampling event. This makes multi-block
#' seasons (e.g. two deployment windows per year with a real gap between
#' them) work correctly without special handling: a day with no recorded
#' event simply has no \code{0} written for it, so it stays \code{NA}
#' everywhere, correctly reflecting "no effort" rather than being forced
#' into a false continuous window.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @param det_df A data frame of detections, one row per (site, species,
#'   date) detected. Must contain the columns named in \code{det_site_col},
#'   \code{det_species_col}, and \code{det_date_col}.
#' @param events_df A data frame of all sampling events (e.g. every
#'   recording), one row per (site, date) of effort — species-agnostic.
#'   Must contain the columns named in \code{events_site_col} and
#'   \code{events_date_col}. Multiple rows per site-day (e.g. several
#'   recordings the same day) are fine; only the presence of effort on
#'   that day matters.
#' @param det_site_col Character. Name of the site column in
#'   \code{det_df}. Defaults to \code{"sites"}.
#' @param det_species_col Character. Name of the species column in
#'   \code{det_df}. Defaults to \code{"species"}.
#' @param det_date_col Character. Name of the date column in
#'   \code{det_df}. Defaults to \code{"date"}.
#' @param events_site_col Character. Name of the site column in
#'   \code{events_df}. Defaults to \code{"site"}.
#' @param events_date_col Character. Name of the date column in
#'   \code{events_df}. Defaults to \code{"eventDate"}.
#'
#' @return A named list with the same shape as
#'   \code{\link{build_detection_array}}, so it is a drop-in replacement:
#' \describe{
#'   \item{y}{3D array \code{[species x sites x dates]}. \code{NA} = no
#'     sampling effort, \code{0} = effort with no detection, \code{1} =
#'     species detected.}
#'   \item{species}{Character vector of species names.}
#'   \item{sites}{Character vector of site names.}
#'   \item{date_seq}{Date vector of the full temporal sequence, spanning
#'     \code{min(events_df$date)} to \code{max(events_df$date)}.}
#' }
#'
#' @details
#' Detections whose site does not appear in \code{events_df}, or whose
#' date falls outside the range covered by \code{events_df}, cannot be
#' placed in the array (there is no corresponding effort cell) and are
#' dropped with a warning rather than causing the function to fail.
#'
#' @seealso \code{\link{build_detection_array}},
#'   \code{\link{build_reindexed_visits}}, \code{\link{build_detection_data}}
#'
#' @importFrom dplyr rename mutate distinct filter
#'
#' @export
build_detection_array_from_events <- function(det_df,
                                              events_df,
                                              det_site_col    = "sites",
                                              det_species_col = "species",
                                              det_date_col    = "date",
                                              events_site_col = "site",
                                              events_date_col = "eventDate") {

  # ── Input validation ─────────────────────────────────────────────────────────
  stopifnot(
    "det_df must be a data frame"    = is.data.frame(det_df),
    "events_df must be a data frame" = is.data.frame(events_df)
  )

  missing_det <- setdiff(c(det_site_col, det_species_col, det_date_col),
                         names(det_df))
  if (length(missing_det) > 0) {
    stop("Missing columns in det_df: ", paste(missing_det, collapse = ", "))
  }

  missing_events <- setdiff(c(events_site_col, events_date_col),
                            names(events_df))
  if (length(missing_events) > 0) {
    stop("Missing columns in events_df: ", paste(missing_events, collapse = ", "))
  }

  # ── Standardise column names and types ─────────────────────────────────────
  det_df <- det_df %>%
    dplyr::rename(
      sites   = !!det_site_col,
      species = !!det_species_col,
      date    = !!det_date_col
    ) %>%
    dplyr::mutate(
      sites   = as.character(sites),
      species = as.character(species),
      date    = as.Date(date)
    )

  events_df <- events_df %>%
    dplyr::rename(
      sites = !!events_site_col,
      date  = !!events_date_col
    ) %>%
    dplyr::mutate(sites = as.character(sites), date = as.Date(date)) %>%
    # one row per site-day of effort, species-agnostic
    dplyr::distinct(sites, date)

  species    <- sort(unique(det_df$species))
  sites      <- sort(unique(events_df$sites))
  full_dates <- seq(min(events_df$date), max(events_df$date), by = "day")

  # ── Drop detections with no corresponding effort cell ──────────────────────
  # A detection whose site never appears in events_df, or whose date falls
  # outside the range covered by events_df, has no valid effort cell to be
  # placed in — keeping it would crash the matrix-index assignment below.
  invalid_site <- !det_df$sites %in% sites
  invalid_date <- det_df$date < min(full_dates) | det_df$date > max(full_dates)
  invalid_rows <- invalid_site | invalid_date

  if (any(invalid_rows)) {
    warning(
      sum(invalid_rows), " detection(s) dropped: site not present in ",
      "events_df or date outside the events_df date range.\n",
      "Affected sites: ",
      paste(unique(det_df$sites[invalid_rows]), collapse = ", ")
    )
    det_df <- det_df[!invalid_rows, ]
  }

  # ── Initialise array ──────────────────────────────────────────────────────────
  y <- array(
    NA,
    dim      = c(length(species), length(sites), length(full_dates)),
    dimnames = list(species = species, sites = sites,
                    date = as.character(full_dates))
  )

  # ── Mark effort: 0 on every (species, site, date) where a recording occurred ──
  idx_site <- match(events_df$sites, sites)
  idx_date <- match(events_df$date, full_dates)

  for (i in seq_len(nrow(events_df))) {
    y[, idx_site[i], idx_date[i]] <- 0
  }

  # ── Overwrite with 1 where a detection actually occurred ───────────────────
  idx_sp     <- match(det_df$species, species)
  idx_site_d <- match(det_df$sites, sites)
  idx_date_d <- match(det_df$date, full_dates)

  y[cbind(idx_sp, idx_site_d, idx_date_d)] <- 1

  list(y = y, species = species, sites = sites, date_seq = full_dates)
}
