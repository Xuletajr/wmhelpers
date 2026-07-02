#' Build a Deployment Table from a Recording Summary Object
#'
#' Constructs a standardised deployment table from a recording summary data
#' frame already loaded in the R environment. The output format is identical
#' to that of \code{\link{read_deployment_file}}, enabling both functions to
#' be used interchangeably in downstream workflows.
#'
#' This function is intended for workflows where deployment data originates
#' from a dashboard export (e.g., \code{event.csv}) rather than a standalone
#' deployment CSV file. The recording summary object is created by grouping
#' events by site and date and counting recordings per site-day.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @param recordings_df A data frame or tibble with one row per site × date
#'   combination. Must contain the columns specified in \code{site_col},
#'   \code{date_col}, \code{lat_col}, \code{lon_col}, and
#'   \code{count_col}. Typically derived from an \code{events} object as:
#'   \preformatted{
#'   recordings_df <- events %>%
#'     dplyr::filter(methods == "audio") %>%
#'     dplyr::group_by(site, eventDate, latitude, longitude) %>%
#'     dplyr::summarise(recording_count = n(), .groups = "drop") %>%
#'     dplyr::arrange(site, eventDate)
#'   }
#' @param method Character. Survey method label used to construct
#'   \code{point_name} and populate the \code{method} column in the output
#'   (e.g., \code{"audio"}, \code{"camera"}, \code{"dna"}).
#' @param site_col Character. Name of the column containing site identifiers.
#'   Site names are converted to lowercase. Defaults to \code{"site"}.
#' @param date_col Character. Name of the column containing sampling dates
#'   (must be \code{Date} class). Defaults to \code{"eventDate"}.
#' @param lat_col Character. Name of the latitude column.
#'   Defaults to \code{"latitude"}.
#' @param lon_col Character. Name of the longitude column.
#'   Defaults to \code{"longitude"}.
#' @param count_col Character. Name of the recording count column. Used
#'   together with \code{min_recordings} to exclude low-effort days before
#'   computing deployment start and end dates. Defaults to
#'   \code{"recording_count"}.
#' @param min_recordings Numeric or \code{NULL}. Minimum number of recordings
#'   required to include a site-day when computing deployment start and end
#'   dates. Site-days below this threshold are excluded before summarising.
#'   Use \code{NULL} (default) to include all days without filtering.
#' @param valid_sites Optional character vector of site names to retain.
#'   Applied after lowercase conversion. Defaults to \code{NULL}.
#' @param issue_sites Optional character vector of site names to exclude.
#'   Applied after lowercase conversion. Defaults to \code{NULL}.
#'
#' @return A tibble with one row per unique site and the following columns,
#'   matching the output of \code{\link{read_deployment_file}}:
#' \describe{
#'   \item{\code{sites}}{Site identifier (lowercase).}
#'   \item{\code{point_name}}{Concatenation of \code{sites} and
#'     \code{method}, separated by \code{"_"}.}
#'   \item{\code{longitude}}{Longitude coordinate (mean across dates if
#'     coordinates vary slightly per event).}
#'   \item{\code{latitude}}{Latitude coordinate (mean across dates).}
#'   \item{\code{start_date}}{Earliest sampling date for the site (after
#'     low-effort filtering if \code{min_recordings} is set).}
#'   \item{\code{end_date}}{Latest sampling date for the site (after
#'     low-effort filtering).}
#'   \item{\code{method}}{Survey method label as supplied in \code{method}.}
#' }
#'
#' @details
#' Site names are normalised to lowercase before any filtering or grouping.
#' Coordinates are averaged per site to handle minor floating-point variation
#' across events at the same location.
#'
#' When \code{min_recordings} is set, site-days with fewer recordings than
#' the threshold are excluded before computing \code{start_date} and
#' \code{end_date}. This ensures that deployment periods reflect only days
#' with sufficient sampling effort, consistent with the behaviour of
#' \code{\link{read_deployment_file}} when \code{min_recordings} is provided.
#'
#' The output is sorted by \code{point_name}, matching
#' \code{\link{read_deployment_file}}.
#'
#' @seealso \code{\link{read_deployment_file}},
#'   \code{\link{set_low_effort_days_to_na}}
#'
#' @importFrom dplyr filter mutate group_by summarise arrange
#' @importFrom rlang sym abort
#'
#' @export
build_deployment_table <- function(recordings_df,
                                   method,
                                   site_col       = "site",
                                   date_col       = "eventDate",
                                   lat_col        = "latitude",
                                   lon_col        = "longitude",
                                   count_col      = "recording_count",
                                   min_recordings = NULL,
                                   valid_sites    = NULL,
                                   issue_sites    = NULL) {

  # ── Input validation ─────────────────────────────────────────────────────────
  stopifnot(
    "recordings_df must be a data frame"   = is.data.frame(recordings_df),
    "method must be a single string"       = is.character(method) &&
      length(method) == 1,
    "site_col must be a single string"     = is.character(site_col) &&
      length(site_col) == 1,
    "date_col must be a single string"     = is.character(date_col) &&
      length(date_col) == 1,
    "lat_col must be a single string"      = is.character(lat_col) &&
      length(lat_col) == 1,
    "lon_col must be a single string"      = is.character(lon_col) &&
      length(lon_col) == 1,
    "count_col must be a single string"    = is.character(count_col) &&
      length(count_col) == 1
  )

  required_cols <- c(site_col, date_col, lat_col, lon_col)
  missing_cols  <- setdiff(required_cols, names(recordings_df))

  if (length(missing_cols) > 0) {
    rlang::abort(paste(
      "The following required columns are missing from recordings_df:",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (!is.null(min_recordings) && !count_col %in% names(recordings_df)) {
    rlang::abort(paste0(
      "Column '", count_col, "' not found in recordings_df. ",
      "It is required when min_recordings is set."
    ))
  }

  # ── Normalise site names to lowercase ─────────────────────────────────────
  df <- recordings_df %>%
    dplyr::mutate(sites = tolower(.data[[site_col]]))

  # ── Apply site filters ────────────────────────────────────────────────────
  if (!is.null(valid_sites)) {
    df <- df %>% dplyr::filter(sites %in% tolower(valid_sites))
  }

  if (!is.null(issue_sites)) {
    df <- df %>% dplyr::filter(!sites %in% tolower(issue_sites))
  }

  # ── Filter low-effort days before computing deployment dates ──────────────
  if (!is.null(min_recordings)) {

    n_before <- nrow(df)

    df <- df %>%
      dplyr::filter(.data[[count_col]] >= min_recordings)

    n_removed <- n_before - nrow(df)

    if (n_removed > 0) {
      message(sprintf(
        "Removed %d site-day(s) with fewer than %d recordings.",
        n_removed, min_recordings
      ))
    }
  }

  # ── Summarise to one row per site ─────────────────────────────────────────
  df %>%
    dplyr::group_by(sites) %>%
    dplyr::summarise(
      longitude  = mean(.data[[lon_col]], na.rm = TRUE),
      latitude   = mean(.data[[lat_col]], na.rm = TRUE),
      start_date = min(.data[[date_col]], na.rm = TRUE),
      end_date   = max(.data[[date_col]], na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    dplyr::mutate(
      method     = method,
      point_name = paste0(sites, "_", method)
    ) %>%
    dplyr::select(
      sites, point_name, longitude, latitude,
      start_date, end_date, method
    ) %>%
    dplyr::arrange(point_name)
}
