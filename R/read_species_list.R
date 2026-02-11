#' Read and standardize species list from a Project Tracking sheet
#'
#' Reads the "Species List" sheet from a Project Tracking Google Sheets
#' document and returns a standardized species metadata table ready for
#' ecological analyses.
#'
#' The function harmonizes column names, removes missing species entries,
#' standardizes species names (first two taxonomic words), and optionally
#' adds a logical marker column.
#'
#' @param sheet_url Character string. URL of the Google Sheets document.
#' @param sheet_name Character string. Name of the sheet to read.
#'   Defaults to `"Species List"`.
#' @param mark Optional character string. If provided, a new column with this
#'   name will be added and filled with `"yes"`.
#'
#' @return A tibble containing a standardized species metadata table.
#'
#' @details
#' Expected columns in the sheet include:
#' \describe{
#'   \item{species}{Scientific name (may include subspecies; only genus and species retained).}
#'   \item{ea}{Flag indicating inclusion in ecological analyses.}
#'   \item{trophic niche}{Species trophic guild classification.}
#'   \item{common name}{English common name.}
#'   \item{habitat density}{Habitat density classification.}
#' }
#'
#' Column names are converted to lower case before processing.
#'
#' @author Thiago Bicudo
#'
#' @export
read_species_list <- function(sheet_url, sheet_name = "Species List", mark = NULL) {
  googlesheets4::read_sheet(sheet_url, sheet = sheet_name) %>%
    dplyr::rename_with(tolower) %>%
    tidyr::drop_na(species) %>%
    dplyr::rename(
      to_analyses   = ea,
      guild         = "trophic niche",
      common_name   = "common name",
      habitatDensity = "habitat density"
    ) %>%
    dplyr::mutate(
      species = stringr::word(species, 1, 2)
    ) %>%
    {
      if (!is.null(mark)) dplyr::mutate(., !!mark := "yes") else .
    } %>%
    dplyr::arrange(species)
}
