#' Coalesce columns with multiple method suffixes into base columns
#'
#' After joining datasets from different monitoring methods (e.g., audio,
#' camera trap, eDNA, citizen science), duplicated columns often appear
#' with suffixes such as `.audio`, `.cam`, `.edna`, etc.
#'
#' This function identifies columns sharing the same base name and
#' coalesces them into a single column, prioritizing non-missing values
#' from left to right according to the order of `suffixes`.
#'
#' @param df A data frame.
#' @param suffixes Character vector of suffixes to detect and merge.
#'   Defaults to `c(".audio", ".cam", ".edna", ".citizen")`.
#'   The order determines coalescing priority.
#'
#' @return A data frame with coalesced base columns and original
#'   suffixed columns removed.
#'
#' @details
#' If `suffixes = c(".audio", ".cam", ".edna")`, the function will:
#' - Identify columns like `species.audio`, `species.cam`, `species.edna`, `species.citizen`
#' - Create a unified column `species`
#' - Use `dplyr::coalesce()` in the order provided
#' - Remove the original suffixed columns
#'
#' @author Thiago Bicudo
#'
#' @export
coalesce_suffix_pairs <- function(df, suffixes = c(".audio", ".cam", ".edna", ".citizen")) {

  pattern <- paste0("(", paste0(suffixes, collapse = "|"), ")$")

  suffixed_cols <- grep(pattern, names(df), value = TRUE)

  if (length(suffixed_cols) == 0) {
    return(df)
  }

  base_names <- unique(sub(pattern, "", suffixed_cols))

  for (base in base_names) {

    candidate_cols <- paste0(base, suffixes)
    existing_cols <- candidate_cols[candidate_cols %in% names(df)]

    if (length(existing_cols) > 0) {
      df[[base]] <- dplyr::coalesce(!!!df[existing_cols])
    }
  }

  dplyr::select(df, -dplyr::matches(pattern))
}
