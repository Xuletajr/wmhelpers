#' Summarise covariate marginal effects from species-level parameter estimates
#'
#' Filters a species-level parameter table by a given covariate prefix and
#' returns the range of credible interval values and a summary of effect
#' directions grouped by one or more trait variables. Designed for outputs of
#' \code{build_param_table()} joined with a species traits table. For
#' integrated multi-method models, pass the data for each detection method
#' separately, as files are exported per method by \code{build_param_table()}.
#'
#' @author Thiago Bicudo
#'
#' @param data A data frame containing species-level parameter estimates,
#'   typically produced by \code{build_param_table()} and joined with a
#'   species traits table via \code{dplyr::left_join()}. Must contain columns
#'   \code{param_spp_code}, \code{CRI_low}, \code{CRI_high}, and \code{Effect},
#'   plus any columns named in \code{group_vars}.
#' @param covariate Character. The covariate prefix to filter by. Must match
#'   the beginning of values in \code{param_spp_code} before the \code{"-"}
#'   separator (e.g., \code{"days"}, \code{"canopyHeight"}, \code{"I(days^2)"}).
#'   Special regex characters are escaped internally, so the string should be
#'   provided as it appears in the data — not as a regex pattern.
#' @param group_vars Character vector. Column names in \code{data} to group
#'   by when summarising effect directions. \code{"Effect"} should always be
#'   included. Defaults to \code{c("guild", "Effect")}. Other useful options
#'   include \code{"habitat"}, \code{"habitatDensity"}, \code{"iucn"}, or any
#'   other trait column available in \code{data}.
#'
#'   Examples:
#'   \preformatted{
#'   group_vars = c("guild", "Effect")           # default
#'   group_vars = c("habitat", "Effect")
#'   group_vars = c("iucn", "Effect")
#'   group_vars = c("guild", "habitat", "Effect")
#'   }
#'
#' @return A named list with three elements: \code{min} (numeric, minimum
#'   value of \code{CRI_low} across all species for the target covariate),
#'   \code{max} (numeric, maximum value of \code{CRI_high}), and
#'   \code{summary} (a \code{tibble} with one row per combination of
#'   \code{group_vars} values, plus a column \code{n} with the species count).
#'
#' @details
#' Filtering uses an anchored regex pattern (\code{"^<escaped_covariate>-"})
#' to match only the exact covariate prefix. Special regex characters in the
#' covariate name (e.g., parentheses and caret in \code{"I(days^2)"}) are
#' automatically escaped before pattern construction, preventing silent
#' mismatches.
#'
#' The anchored pattern also prevents partial matches when covariate names
#' share substrings (e.g., \code{"agb"} will not match
#' \code{"int_agb_ndvi-Species"}).
#'
#' For integrated multi-method models (e.g., audio + camera trap), call this
#' function separately for each method using the corresponding data file
#' exported by \code{build_param_table()}:
#'
#' \preformatted{
#' # Audio
#' summarise_covariate_effects(data = alpha_cov_aud, covariate = "days")
#'
#' # Camera trap
#' summarise_covariate_effects(data = alpha_cov_cam, covariate = "numD")
#' }
#'
#' @importFrom dplyr filter group_by summarise n across all_of
#' @importFrom stringr str_detect str_escape
#'
#' @export
summarise_covariate_effects <- function(data,
                                        covariate,
                                        group_vars = c("guild", "Effect")) {

  # ── Input validation ───────────────────────────────────────────────────────
  stopifnot(
    "data must be a data frame"              = is.data.frame(data),
    "covariate must be a single string"      = is.character(covariate) &&
      length(covariate) == 1 &&
      nchar(covariate) > 0,
    "group_vars must be a character vector"  = is.character(group_vars) &&
      length(group_vars) >= 1
  )

  # ── Check required fixed columns ───────────────────────────────────────────
  required_cols <- c("param_spp_code", "CRI_low", "CRI_high", "Effect")
  missing_cols  <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(
      "The following required columns are missing from data:\n",
      paste(missing_cols, collapse = ", "), "\n",
      "Ensure data was produced by build_param_table() and joined with ",
      "a species traits table."
    )
  }

  # ── Check group_vars columns ───────────────────────────────────────────────
  missing_group <- setdiff(group_vars, names(data))

  if (length(missing_group) > 0) {
    stop(
      "The following group_vars columns are missing from data:\n",
      paste(missing_group, collapse = ", "), "\n",
      "Check that the traits table was joined correctly and contains ",
      "these columns."
    )
  }

  # ── Build anchored regex pattern ───────────────────────────────────────────
  # Special regex characters (e.g. parentheses and ^ in "I(days^2)") are
  # escaped to prevent silent mismatches or errors.
  covariate_pattern <- paste0("^", stringr::str_escape(covariate), "-")

  # ── Filter by covariate ────────────────────────────────────────────────────
  covariate_data <- data %>%
    dplyr::filter(stringr::str_detect(param_spp_code, covariate_pattern))

  # ── Check: covariate found ─────────────────────────────────────────────────
  if (nrow(covariate_data) == 0) {
    available_prefixes <- unique(sub("-.*$", "", data$param_spp_code))
    stop(
      "No rows found in data matching covariate '", covariate, "'.\n",
      "Expected pattern: '", covariate_pattern,
      "' in column 'param_spp_code'.\n",
      "Available prefixes: ",
      paste(available_prefixes, collapse = ", ")
    )
  }

  # ── Check: NA in CRI columns ───────────────────────────────────────────────
  n_na_low  <- sum(is.na(covariate_data$CRI_low))
  n_na_high <- sum(is.na(covariate_data$CRI_high))

  if (n_na_low > 0 || n_na_high > 0) {
    warning(
      n_na_low, " NA(s) in CRI_low and ",
      n_na_high, " NA(s) in CRI_high for covariate '", covariate, "'. ",
      "These will be ignored in min/max computation."
    )
  }

  # ── Check: NA in group_vars columns ───────────────────────────────────────
  for (gv in group_vars) {
    n_na <- sum(is.na(covariate_data[[gv]]))
    if (n_na > 0) {
      warning(
        n_na, " species with NA in '", gv, "' for covariate '", covariate,
        "'. These will appear as NA in the summary table. ",
        "Check that the traits join was complete."
      )
    }
  }

  # ── Compute CRI range ──────────────────────────────────────────────────────
  min_val <- min(covariate_data$CRI_low,  na.rm = TRUE)
  max_val <- max(covariate_data$CRI_high, na.rm = TRUE)

  # ── Summarise by group_vars and effect direction ───────────────────────────
  summary_data <- covariate_data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  # ── Check: summary not empty ───────────────────────────────────────────────
  if (nrow(summary_data) == 0) {
    warning(
      "Summary table is empty for covariate '", covariate, "'. ",
      "Check that group_vars columns contain valid values."
    )
  }

  return(list(
    min     = min_val,
    max     = max_val,
    summary = summary_data
  ))
}
