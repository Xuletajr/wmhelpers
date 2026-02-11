#' Read and Standardise Deployment Data (Camera or Audio)
#'
#' Reads deployment data from camera trap or audio surveys and returns
#' a standardised data frame with consistent structure across methods.
#'
#' For camera data, start and end dates must already exist in the file.
#' For audio data, year/month/day columns are combined to compute
#' start and end dates per site, and coordinates are joined externally.
#'
#' @param path Character. Path to deployment CSV file.
#' @param site_col Character. Column name containing site identifiers.
#' @param method Character. Either `"camera"` or `"audio"`.
#' @param valid_sites Optional character vector of valid site names to retain.
#' @param issue_sites Optional character vector of site names to exclude.
#' @param coord_path Optional character. Path to coordinate CSV file.
#' Required when `method = "audio"`.
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
#' @examples
#' \dontrun{
#' read_deployment_file(
#'   path = "deployments.csv",
#'   site_col = "deployment_id",
#'   method = "camera"
#' )
#' }
#'
#'#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @export
read_deployment_file <- function(path,
                                 site_col,
                                 method = c("camera", "audio"),
                                 valid_sites = NULL,
                                 issue_sites = NULL,
                                 coord_path = NULL) {

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

    coords <- readr::read_csv(coord_path,
                              show_col_types = FALSE)

    df <- df %>%
      dplyr::mutate(
        date = lubridate::make_date(year, month, day)
      ) %>%
      dplyr::rename(sites = !!rlang::sym(site_col)) %>%
      dplyr::left_join(coords, by = c("sites" = "site_name")) %>%
      dplyr::group_by(sites, longitude, latitude) %>%
      dplyr::summarise(
        start_date = min(date),
        end_date   = max(date),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        sites = tolower(sites)
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
