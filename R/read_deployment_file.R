#' Read and Standardise Deployment Data (Camera or Audio)
#'
#' Reads deployment data from camera trap or audio surveys and returns
#' a standardised data frame with consistent structure across methods.
#'
#' For camera data, start and end dates must already exist in the file.
#' For audio data, year/month/day columns are combined to compute
#' start and end dates per site. Days with fewer recordings than
#' `min_recordings` are excluded before calculating deployment periods.
#'
#' @param path Character. Path to deployment CSV file.
#' @param site_col Character. Column name containing site identifiers.
#' @param method Character. Either `"camera"` or `"audio"`.
#' @param valid_sites Optional character vector of valid site names to retain.
#' @param issue_sites Optional character vector of site names to exclude.
#' @param coord_path Optional character. Path to coordinate CSV file.
#'   Required when `method = "audio"`.
#' @param min_recordings Optional numeric. Minimum number of recordings
#'   required to keep a site-day. Days below this threshold are excluded
#'   before computing start/end dates. Only used when `method = "audio"`.
#'   Requires a `recording_count` column in the input file. Default is NULL
#'   (no filtering).
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{sites}{Site identifier (lowercase)}
#'   \item{point_name}{Concatenation of site and method}
#'   \item{longitude}{Longitude coordinate}
#'   \item{latitude}{Latitude coordinate}
#'   \item{start_date}{Deployment start date}
#'   \item{end_date}{Deployment end date}
#'   \item{method}{Survey method}
#' }
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @export
read_deployment_file3 <- function(path,
                                  site_col,
                                  method = c("camera", "audio"),
                                  valid_sites = NULL,
                                  issue_sites = NULL,
                                  coord_path = NULL,
                                  min_recordings = NULL) {

  method <- match.arg(method)

  if (!file.exists(path)) {
    rlang::abort(paste("File not found:", path))
  }

  df <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE),
    error = function(e) {
      rlang::abort(paste("Failed to read:", path))
    }
  )

  # -------------------------
  # CAMERA METHOD
  # -------------------------
  if (method == "camera") {

    required_cols <- c(site_col, "longitude",
                       "latitude", "start_date", "end_date")

    missing_cols <- setdiff(required_cols, names(df))

    if (length(missing_cols) > 0) {
      rlang::abort(paste("Missing columns:",
                         paste(missing_cols, collapse = ", ")))
    }

    df <- df %>%
      dplyr::rename(sites = !!rlang::sym(site_col)) %>%
      dplyr::mutate(
        sites      = tolower(sites),
        start_date = as.Date(start_date),
        end_date   = as.Date(end_date)
      ) %>%
      dplyr::select(
        sites, longitude, latitude,
        start_date, end_date
      )
  }

  # -------------------------
  # AUDIO METHOD
  # -------------------------
  if (method == "audio") {

    required_cols <- c(site_col, "year", "month", "day")
    missing_cols <- setdiff(required_cols, names(df))

    if (length(missing_cols) > 0) {
      rlang::abort(paste("Missing columns:",
                         paste(missing_cols, collapse = ", ")))
    }

    if (is.null(coord_path) || !file.exists(coord_path)) {
      rlang::abort("coord_path must be provided for audio method.")
    }

    # Normalise site names to lowercase BEFORE joining
    coords <- readr::read_csv(coord_path, show_col_types = FALSE) %>%
      dplyr::mutate(site_name = tolower(site_name))

    df <- df %>%
      dplyr::rename(sites = !!rlang::sym(site_col)) %>%
      dplyr::mutate(
        sites = tolower(sites),
        date  = lubridate::make_date(year, month, day)
      )

    # ------------------------------------
    # FILTER LOW-EFFORT DAYS (by threshold)
    # ------------------------------------
    if (!is.null(min_recordings)) {

      if (!"recording_count" %in% names(df)) {
        rlang::abort(
          "Column 'recording_count' not found. It is required when min_recordings is set."
        )
      }

      n_before <- nrow(df)

      df <- df %>%
        dplyr::filter(recording_count >= min_recordings)

      n_removed <- n_before - nrow(df)
      message(sprintf("Removed %d site-day(s) with fewer than %d recordings.",
                      n_removed, min_recordings))
    }

    df <- df %>%
      dplyr::left_join(coords, by = c("sites" = "site_name")) %>%
      dplyr::group_by(sites, longitude, latitude) %>%
      dplyr::summarise(
        start_date = min(date),
        end_date   = max(date),
        .groups    = "drop"
      )
  }

  # -------------------------
  # COMMON PROCESSING
  # -------------------------
  df <- df %>%
    dplyr::mutate(
      method     = method,
      point_name = paste0(sites, "_", method)
    )

  if (!is.null(valid_sites)) {
    df <- df %>% dplyr::filter(sites %in% valid_sites)
  }

  if (!is.null(issue_sites)) {
    df <- df %>% dplyr::filter(!sites %in% issue_sites)
  }

  df %>%
    dplyr::arrange(point_name)
}
