#' Summarise covariate marginal effects from species-level parameter estimates
#'
#' Filters a species-level parameter table by a given covariate prefix and
#' returns the range of credible interval values and a summary of effect
#' directions by guild. Designed for outputs of \code{build_param_table()}
#' joined with a species traits table. For integrated multi-method models,
#' pass the data for each detection method separately, as files are exported
#' per method by \code{build_param_table()}.
#'
#' @author Thiago Bicudo
#'
#' @param data A data frame containing species-level parameter estimates,
#'   typically produced by \code{build_param_table()} and joined with a
#'   species traits table via \code{dplyr::left_join()}. Must contain columns
#'   \code{param_spp_code}, \code{CRI_low}, \code{CRI_high}, \code{guild},
#'   and \code{Effect}.
#' @param covariate Character. The covariate prefix to filter by. Must match
#'   the beginning of values in \code{param_spp_code} before the \code{"-"}
#'   separator (e.g., \code{"days"}, \code{"canopyHeight"}, \code{"I(days^2)"}).
#'   Special regex characters are escaped internally, so the string should be
#'   provided as it appears in the data — not as a regex pattern.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{\code{min}}{Numeric. Minimum value of \code{CRI_low} across all
#'     species for the target covariate.}
#'   \item{\code{max}}{Numeric. Maximum value of \code{CRI_high} across all
#'     species for the target covariate.}
#'   \item{\code{summary}}{A \code{tibble} with columns \code{guild},
#'     \code{Effect}, and \code{n}, summarising the number of species per
#'     guild and effect direction.}
#' }
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
#' @importFrom dplyr filter group_by summarise n
#' @importFrom stringr str_detect str_escape
#'
#' @export
summarise_covariate_effects <- function(data, covariate) {

  # ── Input validation ───────────────────────────────────────────────────────
  stopifnot(
    "data must be a data frame"         = is.data.frame(data),
    "covariate must be a single string" = is.character(covariate) &&
      length(covariate) == 1 &&
      nchar(covariate) > 0
  )

  required_cols <- c("param_spp_code", "CRI_low", "CRI_high", "guild", "Effect")
  missing_cols  <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(
      "The following required columns are missing from data:\n",
      paste(missing_cols, collapse = ", "), "\n",
      "Ensure data was produced by build_param_table() and joined with ",
      "a species traits table containing 'guild'."
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
    available_prefixes <- unique(
      sub("-.*$", "", data$param_spp_code)
    )
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
      n_na_low,  " NA(s) in CRI_low and ",
      n_na_high, " NA(s) in CRI_high for covariate '", covariate, "'. ",
      "These will be ignored in min/max computation."
    )
  }

  # ── Check: NA in guild column ──────────────────────────────────────────────
  n_na_guild <- sum(is.na(covariate_data$guild))

  if (n_na_guild > 0) {
    warning(
      n_na_guild, " species with NA in 'guild' for covariate '", covariate,
      "'. These will appear as NA in the summary table. ",
      "Check that the traits join was complete."
    )
  }

  # ── Compute CRI range ──────────────────────────────────────────────────────
  min_val <- min(covariate_data$CRI_low,  na.rm = TRUE)
  max_val <- max(covariate_data$CRI_high, na.rm = TRUE)

  # ── Summarise by guild and effect direction ────────────────────────────────
  summary_data <- covariate_data %>%
    dplyr::group_by(guild, Effect) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  # ── Check: summary not empty ───────────────────────────────────────────────
  if (nrow(summary_data) == 0) {
    warning(
      "Summary table is empty for covariate '", covariate, "'. ",
      "Check that 'guild' and 'Effect' columns contain valid values."
    )
  }

  return(list(
    min     = min_val,
    max     = max_val,
    summary = summary_data
  ))
}
