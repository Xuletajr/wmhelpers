#' Update Appendix with Occupancy and Detection Probabilities
#'
#' Joins processed occupancy and detection probability estimates to an existing
#' appendix data frame. Supports two modes: a simple intercept-only model
#' (single probability column per species) or a class/year model (one column
#' per class level).
#'
#' @param appendix A data frame containing at minimum a column named
#'   \code{Species}, to which occupancy and detection columns will be joined.
#' @param occ_data A data frame of occupancy estimates. Must contain columns
#'   \code{species} and \code{mean}, and — when \code{use_classes = TRUE} —
#'   a column \code{class}.
#' @param det_data A data frame of detection estimates. Must contain columns
#'   \code{species} and \code{mean}, and — when \code{use_classes = TRUE} —
#'   a column \code{class}.
#' @param use_classes Logical. If \code{TRUE}, probabilities are pivoted into
#'   one column per class level. Defaults to \code{FALSE}.
#' @param occ_pattern Optional character string. Overrides default naming
#'   convention for occupancy columns when \code{use_classes = TRUE}.
#' @param det_pattern Optional character string. Overrides default naming
#'   convention for detection columns when \code{use_classes = TRUE}.
#'
#' @return A data frame equivalent to \code{appendix} with additional columns
#'   for occupancy and detection probabilities (rounded to 3 decimal places).
#'   The \code{Species} column is always returned first.
#'
#' @details
#' When \code{use_classes = FALSE}, the output columns are named
#' \code{Occupancy.probability} and \code{Detection.probability}.
#'
#' When \code{use_classes = TRUE}, column names are derived from the
#' \code{class} variable (e.g. \code{occ_prob_2020}, \code{det_prob_forest}),
#' unless overridden by \code{occ_pattern} or \code{det_pattern}.
#'
#' @importFrom dplyr ungroup mutate select rename left_join everything across
#' @importFrom tidyr pivot_wider
#'
#' @examples
#' \dontrun{
#' # Intercept-only model
#' updated <- update_appendix(
#'   appendix = my_appendix,
#'   occ_data = occ_estimates,
#'   det_data = det_estimates
#' )
#'
#' # Class-level model
#' updated <- update_appendix(
#'   appendix    = my_appendix,
#'   occ_data    = occ_estimates,
#'   det_data    = det_estimates,
#'   use_classes = TRUE
#' )
#' }
#'
#' @author Thiago Bicudo
#' @export
update_appendix <- function(appendix,
                            occ_data,
                            det_data,
                            use_classes = FALSE,
                            occ_pattern = NULL,
                            det_pattern = NULL) {

  # ── Input validation ────────────────────────────────────────────────
  stopifnot(
    is.data.frame(appendix),
    is.data.frame(occ_data),
    is.data.frame(det_data),
    "Species" %in% names(appendix),
    "species" %in% names(occ_data),
    "species" %in% names(det_data),
    "mean" %in% names(occ_data),
    "mean" %in% names(det_data)
  )

  if (use_classes) {
    stopifnot(
      "class" %in% names(occ_data),
      "class" %in% names(det_data)
    )
  }

  # ── Helper: column naming ───────────────────────────────────────────
  make_col_name <- function(class_vec, prefix, pattern) {
    if (!is.null(pattern)) {
      pattern
    } else {
      paste0(prefix, class_vec)
    }
  }

  # ── Helper: process class data ──────────────────────────────────────
  process_class_data <- function(data, prefix, pattern) {

    data %>%
      ungroup() %>%
      mutate(col_name = make_col_name(class, prefix, pattern)) %>%
      select(species, col_name, mean) %>%
      pivot_wider(
        names_from  = col_name,
        values_from = mean
      ) %>%
      mutate(across(starts_with(prefix), ~ round(.x, 3)))
  }

  # ── Main logic ──────────────────────────────────────────────────────
  if (use_classes) {

    occ_processed <- process_class_data(occ_data, "occ_prob_", occ_pattern)
    det_processed <- process_class_data(det_data, "det_prob_", det_pattern)

  } else {

    occ_processed <- occ_data %>%
      ungroup() %>%
      select(species, mean) %>%
      rename(Occupancy.probability = mean) %>%
      mutate(Occupancy.probability = round(Occupancy.probability, 3))

    det_processed <- det_data %>%
      ungroup() %>%
      select(species, mean) %>%
      rename(Detection.probability = mean) %>%
      mutate(Detection.probability = round(Detection.probability, 3))
  }

  # ── Join ────────────────────────────────────────────────────────────
  appendix %>%
    left_join(occ_processed, by = c("Species" = "species")) %>%
    left_join(det_processed, by = c("Species" = "species")) %>%
    select(Species, everything())
}
