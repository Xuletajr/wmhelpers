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
#' @param use_classes Logical. If \code{TRUE}, the function pivots probabilities
#'   into one column per class level (e.g. \code{occ_prob_forest},
#'   \code{det_prob_urban}). If \code{FALSE} (default), a single
#'   \code{Occupancy.probability} and \code{Detection.probability} column are
#'   added.
#' @param occ_pattern Optional character string. When \code{use_classes = TRUE},
#'   overrides the default column naming convention (\code{paste0("occ_prob_",
#'   class)}) with a fixed column name. Useful when all rows share the same
#'   class label. Default is \code{NULL}.
#' @param det_pattern Optional character string. When \code{use_classes = TRUE},
#'   overrides the default column naming convention (\code{paste0("det_prob_",
#'   class)}) with a fixed column name. Default is \code{NULL}.
#'
#' @return A data frame equivalent to \code{appendix} with additional columns
#'   for occupancy and detection probabilities (rounded to 3 decimal places),
#'   joined by \code{Species}. The \code{Species} column is always returned
#'   first.
#'
#' @details
#' When \code{use_classes = FALSE}, the output columns are named
#' \code{Occupancy.probability} and \code{Detection.probability}.
#'
#' When \code{use_classes = TRUE}, column names are derived from the
#' \code{class} variable in the input data (e.g. \code{occ_prob_2020},
#' \code{det_prob_forest}), unless overridden by \code{occ_pattern} or
#' \code{det_pattern}.
#'
#' Joining is performed via a left join on \code{Species} (in \code{appendix})
#' and \code{species} (in \code{occ_data}/\code{det_data}), so species present
#' in \code{appendix} but absent from the probability data will receive
#' \code{NA} values.
#'
#' @importFrom dplyr ungroup mutate select rename left_join everything
#' @importFrom tidyr pivot_wider
#'
#' @examples
#' \dontrun{
#' # Intercept-only model
#' updated <- update.appendix(
#'   appendix = my_appendix,
#'   occ_data = occ_estimates,
#'   det_data = det_estimates
#' )
#'
#' # Class-level model
#' updated <- update.appendix(
#'   appendix    = my_appendix,
#'   occ_data    = occ_estimates,
#'   det_data    = det_estimates,
#'   use_classes = TRUE
#' )
#' }
#'
#' @export
update.appendix <- function(appendix,
                            occ_data,
                            det_data,
                            use_classes = FALSE,
                            occ_pattern = NULL,
                            det_pattern = NULL) {

  # ── Input validation ────────────────────────────────────────────────────────

  stopifnot(
    "appendix must be a data frame"       = is.data.frame(appendix),
    "occ_data must be a data frame"       = is.data.frame(occ_data),
    "det_data must be a data frame"       = is.data.frame(det_data),
    "'Species' column missing in appendix" = "Species" %in% names(appendix),
    "'species' column missing in occ_data" = "species" %in% names(occ_data),
    "'species' column missing in det_data" = "species" %in% names(det_data),
    "'mean' column missing in occ_data"    = "mean"    %in% names(occ_data),
    "'mean' column missing in det_data"    = "mean"    %in% names(det_data)
  )

  if (use_classes) {
    stopifnot(
      "'class' column required in occ_data when use_classes = TRUE" =
        "class" %in% names(occ_data),
      "'class' column required in det_data when use_classes = TRUE" =
        "class" %in% names(det_data)
    )
  }

  # ── Helper: build column name ────────────────────────────────────────────────
  # Returns either the user-supplied override pattern or the default
  # constructed by pasting a prefix with the class label.

  make_col_name <- function(data, prefix, pattern) {
    if (!is.null(pattern)) {
      pattern                        # fixed override (e.g. "occ_prob_forest")
    } else {
      paste0(prefix, data$class)     # dynamic: "occ_prob_<class>"
    }
  }

  # ── Helper: process a probability data frame ─────────────────────────────────
  # Shared pipeline for both occ_data and det_data when use_classes = TRUE.
  # Attaches a col_name column, keeps only needed columns, widens to one
  # column per class, and rounds all probability columns to 3 d.p.

  process_class_data <- function(data, prefix, pattern) {
    col_prefix <- prefix   # e.g. "occ_prob_" or "det_prob_"

    data %>%
      ungroup() %>%
      mutate(col_name = make_col_name(data, col_prefix, pattern)) %>%
      select(species, col_name, mean) %>%
      pivot_wider(
        names_from  = col_name,
        values_from = mean
      ) %>%
      # Round all newly created probability columns
      mutate(across(starts_with(col_prefix), ~ round(.x, 3)))
  }

  # ── Branch: class-level vs. intercept-only ───────────────────────────────────

  if (use_classes) {

    # One column per class level for both occupancy and detection
    occ_processed <- process_class_data(occ_data, "occ_prob_", occ_pattern)
    det_processed <- process_class_data(det_data, "det_prob_", det_pattern)

  } else {

    # Single occupancy probability column per species
    occ_processed <- occ_data %>%
      ungroup() %>%
      select(species, mean) %>%
      rename(Occupancy.probability = mean) %>%
      mutate(Occupancy.probability = round(Occupancy.probability, 3))

    # Single detection probability column per species
    det_processed <- det_data %>%
      ungroup() %>%
      select(species, mean) %>%
      rename(Detection.probability = mean) %>%
      mutate(Detection.probability = round(Detection.probability, 3))
  }

  # ── Join and return ──────────────────────────────────────────────────────────
  # Left joins preserve all rows in appendix; species absent from the
  # probability data will receive NA in the new columns.

  appendix %>%
    left_join(occ_processed, by = c("Species" = "species")) %>%
    left_join(det_processed, by = c("Species" = "species")) %>%
    select(Species, everything())   # ensure Species stays the first column
}