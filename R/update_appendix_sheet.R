#' Update Appendix with Occupancy and Detection Probabilities
#'
#' Joins processed occupancy and detection probability estimates to an existing
#' appendix data frame. Supports two modes: a simple intercept-only model
#' (single probability column per species) or a class/year model (one column
#' per class level). For integrated multi-method models, where occupancy is
#' shared but detection is estimated separately per method, use
#' \code{skip_occ = TRUE} on subsequent calls to avoid duplicate occupancy
#' columns.
#'
#' @param appendix A data frame containing at minimum a column named
#'   \code{Species}, to which occupancy and detection columns will be joined.
#' @param occ_data A data frame of occupancy estimates. Must contain columns
#'   \code{species} and \code{mean}, and — when \code{use_classes = TRUE} —
#'   a column \code{class}. Ignored when \code{skip_occ = TRUE}.
#' @param det_data A data frame of detection estimates. Must contain columns
#'   \code{species} and \code{mean}, and — when \code{use_classes = TRUE} —
#'   a column \code{class}.
#' @param use_classes Logical. If \code{TRUE}, probabilities are pivoted into
#'   one column per class level. Defaults to \code{FALSE}.
#' @param skip_occ Logical. If \code{TRUE}, the occupancy join is skipped
#'   entirely. Intended for integrated models where \code{update_appendix()}
#'   is called more than once and occupancy has already been added in a
#'   previous call. Defaults to \code{FALSE}.
#' @param occ_pattern Optional character string. Overrides default naming
#'   convention for occupancy columns when \code{use_classes = TRUE}.
#' @param det_pattern Optional character string. Overrides default naming
#'   convention for detection columns when \code{use_classes = TRUE}.
#'
#' @return A data frame equivalent to \code{appendix} with additional columns
#'   for occupancy (unless \code{skip_occ = TRUE}) and detection probabilities,
#'   rounded to 3 decimal places. The \code{Species} column is always returned
#'   first.
#'
#' @details
#' When \code{use_classes = FALSE}, the output columns are named
#' \code{Occupancy.probability} and \code{Detection.probability}.
#'
#' When \code{use_classes = TRUE}, column names are derived from the
#' \code{class} variable (e.g. \code{occ_prob_2020}, \code{det_prob_forest}),
#' unless overridden by \code{occ_pattern} or \code{det_pattern}.
#'
#' For integrated models with multiple detection methods, the recommended
#' workflow is:
#' \preformatted{
#' # First call: occupancy + first detection method
#' appendix2 <- update_appendix(
#'   appendix,
#'   occ_data    = export_occ_prob,
#'   det_data    = export_det_prob_aud,
#'   use_classes = TRUE,
#'   occ_pattern = "occupancy",
#'   det_pattern = "detection_audio"
#' )
#'
#' # Subsequent calls: detection only (skip_occ = TRUE)
#' appendix3 <- update_appendix(
#'   appendix2,
#'   occ_data    = NULL,
#'   det_data    = export_det_prob_cam,
#'   use_classes = TRUE,
#'   skip_occ    = TRUE,
#'   det_pattern = "detection_camera"
#' )
#' }
#'
#' @importFrom dplyr ungroup mutate select rename left_join everything across
#' @importFrom tidyr pivot_wider
#'
#' @author Thiago Bicudo
#' @export
update_appendix_sheet <- function(appendix,
                                  occ_data    = NULL,
                                  det_data,
                                  use_classes = FALSE,
                                  skip_occ    = FALSE,
                                  occ_pattern = NULL,
                                  det_pattern = NULL) {

  # ── Input validation ───────────────────────────────────────────────────────
  stopifnot(
    "appendix must be a data frame"              = is.data.frame(appendix),
    "appendix must contain a 'Species' column"   = "Species" %in% names(appendix),
    "det_data must be a data frame"              = is.data.frame(det_data),
    "det_data must contain a 'species' column"   = "species" %in% names(det_data),
    "det_data must contain a 'mean' column"      = "mean" %in% names(det_data)
  )

  if (!skip_occ) {
    stopifnot(
      "occ_data must be a data frame when skip_occ = FALSE" =
        is.data.frame(occ_data),
      "occ_data must contain a 'species' column" = "species" %in% names(occ_data),
      "occ_data must contain a 'mean' column"    = "mean" %in% names(occ_data)
    )
  }

  if (use_classes) {
    stopifnot(
      "det_data must contain a 'class' column" = "class" %in% names(det_data)
    )
    if (!skip_occ) {
      stopifnot(
        "occ_data must contain a 'class' column when use_classes = TRUE" =
          "class" %in% names(occ_data)
      )
    }
  }

  # ── Helper: column naming ──────────────────────────────────────────────────
  make_col_name <- function(class_vec, prefix, pattern) {
    if (!is.null(pattern)) pattern else paste0(prefix, class_vec)
  }

  # ── Helper: process class data ─────────────────────────────────────────────
  process_class_data <- function(data, prefix, pattern) {
    data %>%
      dplyr::ungroup() %>%
      dplyr::mutate(col_name = make_col_name(class, prefix, pattern)) %>%
      dplyr::select(species, col_name, mean) %>%
      tidyr::pivot_wider(
        names_from  = col_name,
        values_from = mean
      ) %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))
  }

  # ── Process detection ──────────────────────────────────────────────────────
  det_processed <- if (use_classes) {
    process_class_data(det_data, "det_prob_", det_pattern)
  } else {
    det_data %>%
      dplyr::ungroup() %>%
      dplyr::select(species, mean) %>%
      dplyr::rename(Detection.probability = mean) %>%
      dplyr::mutate(Detection.probability = round(Detection.probability, 3))
  }

  # ── Process occupancy (unless skip_occ = TRUE) ────────────────────────────
  occ_processed <- if (!skip_occ) {
    if (use_classes) {
      process_class_data(occ_data, "occ_prob_", occ_pattern)
    } else {
      occ_data %>%
        dplyr::ungroup() %>%
        dplyr::select(species, mean) %>%
        dplyr::rename(Occupancy.probability = mean) %>%
        dplyr::mutate(Occupancy.probability = round(Occupancy.probability, 3))
    }
  }

  # ── Join ───────────────────────────────────────────────────────────────────
  result <- appendix %>%
    dplyr::left_join(det_processed, by = c("Species" = "species"))

  if (!skip_occ) {
    result <- result %>%
      dplyr::left_join(occ_processed, by = c("Species" = "species"))
  }

  result %>% dplyr::select(Species, dplyr::everything())
}
