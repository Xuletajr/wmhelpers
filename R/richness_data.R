# =============================================================================
# Richness extraction and data preparation functions
# File: R/richness_data.R
# =============================================================================


#' Extract species richness estimates from the Z matrix
#'
#' Computes site-level species richness estimates from the posterior samples
#' of the latent occupancy matrix (Z). Each column of the resulting matrix
#' represents the number of species present at a given site for one MCMC sample.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta") and Thiago Bicudo
#'
#' @param z_samples A 3-dimensional array with dimensions
#'   \code{[MCMC x species x sites]}, typically \code{out_bestmodel$z.samples}
#'   from a fitted \code{spOccupancy} model.
#' @param site_names Optional character vector of site names to assign as row
#'   names to the resulting richness matrix. Must match the number of sites
#'   (third dimension of \code{z_samples}). Defaults to \code{NULL}.
#' @param species_idx Optional integer vector of species indices to subset
#'   before computing richness. If \code{NULL} (default), all species are used.
#'   Useful for computing guild- or group-specific richness. Use
#'   \code{get_species_idx()} to derive indices from a traits table.
#'
#' @return A numeric matrix of dimensions \code{[sites x MCMC]} containing
#'   species richness estimates per site per posterior sample.
#'
#' @details
#' Richness at each site is computed as the column sum of the species x sites
#' slice for each MCMC sample. If \code{species_idx} is provided, only the
#' selected species contribute to the richness count — enabling guild- or
#' group-specific analyses.
#'
#' @seealso \code{\link{get_species_idx}}, \code{\link{build_richness_summary}}
#'
#' @export
extract_richness <- function(z_samples,
                             site_names  = NULL,
                             species_idx = NULL) {

  stopifnot(
    "z_samples must be a 3D array"        = is.array(z_samples) &&
      length(dim(z_samples)) == 3,
    "site_names must be character or NULL" = is.null(site_names) ||
      is.character(site_names),
    "species_idx must be integer or NULL"  = is.null(species_idx) ||
      is.numeric(species_idx)
  )

  if (!is.null(species_idx)) {
    z_samples <- z_samples[, species_idx, , drop = FALSE]
  }

  n_mcmc  <- dim(z_samples)[1]
  n_sites <- dim(z_samples)[3]

  if (!is.null(site_names) && length(site_names) != n_sites) {
    stop(
      "Length of site_names (", length(site_names), ") must match the ",
      "number of sites in z_samples (", n_sites, ")."
    )
  }

  richness_matrix <- matrix(NA_real_, nrow = n_sites, ncol = n_mcmc)

  for (i in seq_len(n_mcmc)) {
    richness_matrix[, i] <- colSums(z_samples[i, , ])
  }

  if (!is.null(site_names)) {
    rownames(richness_matrix) <- site_names
  }

  richness_matrix
}


#' Get species indices from a traits table
#'
#' Translates a filter applied to a species traits table into a vector of
#' integer indices aligned with the model species order (\code{sp.names}).
#' The resulting vector can be passed directly to
#' \code{extract_richness(species_idx = ...)}.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta") and Thiago Bicudo
#'
#' @param sp_names Character vector of species names in model order, typically
#'   \code{out_bestmodel$sp.names}.
#' @param traits_df A data frame with at least a species name column and one or
#'   more trait columns to filter on.
#' @param sp_col Character. Name of the species name column in
#'   \code{traits_df}. Defaults to \code{"species"}.
#' @param ... Filtering expressions passed to \code{dplyr::filter()}.
#'   For example, \code{guild == "Frugivore"} or
#'   \code{habitat == "Forest", group == "Aves"}.
#'
#' @return An integer vector of species indices relative to \code{sp_names}.
#'   Species present in the filtered \code{traits_df} but absent from
#'   \code{sp_names} are dropped with a warning. Stops with an error if no
#'   species remain after filtering.
#'
#' @details
#' The function uses \code{...} to pass filtering expressions directly to
#' \code{dplyr::filter()}, allowing any combination of trait conditions.
#' Only species that appear in both the filtered traits table and
#' \code{sp_names} are returned.
#'
#' @examples
#' \dontrun{
#' frugivore_idx <- get_species_idx(
#'   sp_names  = out_bestmodel$sp.names,
#'   traits_df = species_traits,
#'   sp_col    = "species",
#'   guild == "Frugivore"
#' )
#'
#' richness_frugivores <- extract_richness(
#'   z_samples   = out_bestmodel$z.samples,
#'   site_names  = occ.covs$sites,
#'   species_idx = frugivore_idx
#' )
#' }
#'
#' @seealso \code{\link{extract_richness}}
#'
#' @export
get_species_idx <- function(sp_names,
                            traits_df,
                            sp_col = "species",
                            ...) {

  stopifnot(
    "sp_names must be a character vector"  = is.character(sp_names),
    "traits_df must be a data frame"       = is.data.frame(traits_df),
    "sp_col must be a single string"       = is.character(sp_col) &&
      length(sp_col) == 1,
    "sp_col must exist in traits_df"       = sp_col %in% names(traits_df)
  )

  filtered       <- traits_df %>% dplyr::filter(...)
  target_species <- filtered[[sp_col]]

  not_found <- setdiff(target_species, sp_names)
  if (length(not_found) > 0) {
    warning(
      length(not_found), " species in traits_df not found in sp_names ",
      "and will be ignored:\n",
      paste(not_found, collapse = ", ")
    )
  }

  idx <- which(sp_names %in% target_species)

  if (length(idx) == 0) {
    stop(
      "No species matched after filtering. ",
      "Check filter expression and sp_col."
    )
  }

  message(length(idx), " species selected out of ",
          length(sp_names), " in the model.")

  idx
}


#' Build a site-level richness summary table
#'
#' Merges site-level richness estimates (from the Z-matrix posterior samples)
#' with associated site metadata and computes posterior summary statistics
#' (mean, SD, median, and credible intervals) per site and optional grouping
#' variables.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta") and Thiago Bicudo
#'
#' @param richness_matrix A \code{[sites x MCMC]} numeric matrix, output from
#'   \code{extract_richness()}.
#' @param site_data A data frame containing site-level metadata. Must include
#'   a column matching \code{site_col}.
#' @param group_vars Character vector of grouping variable names in
#'   \code{site_data} to include in the output. If \code{NULL}, only site-level
#'   summaries are returned.
#' @param site_col Character. Name of the site identifier column in
#'   \code{site_data}. Must match \code{rownames(richness_matrix)}.
#'   Defaults to \code{"sites"}.
#'
#' @return A \code{tibble} with one row per site, containing the site
#'   identifier, any requested grouping variables, and posterior summary
#'   columns: \code{mean}, \code{sd}, \code{median}, \code{cil} (2.5th
#'   percentile), \code{ciu} (97.5th percentile), \code{q05}, and \code{q95}.
#'
#' @details
#' Row names of \code{richness_matrix} must match values in the
#' \code{site_col} column of \code{site_data}. An informative error is raised
#' if any row names are not found.
#'
#' MCMC sample columns are identified automatically regardless of their names
#' (using \code{tidyr::pivot_longer(cols = -all_of(site_col))}) — the function
#' is robust to different column naming conventions from \code{as_tibble()}.
#'
#' @seealso \code{\link{extract_richness}}, \code{\link{build_richness_density}}
#'
#' @export
build_richness_summary <- function(richness_matrix,
                                   site_data,
                                   group_vars = NULL,
                                   site_col   = "sites") {

  stopifnot(
    "richness_matrix must be a numeric matrix" = is.matrix(richness_matrix),
    "site_data must be a data frame"           = is.data.frame(site_data),
    "site_col must be a single string"         = is.character(site_col) &&
      length(site_col) == 1,
    "site_col must exist in site_data"         = site_col %in% names(site_data)
  )

  missing_sites <- setdiff(rownames(richness_matrix), site_data[[site_col]])
  if (length(missing_sites) > 0) {
    stop(
      "The following rownames from richness_matrix were not found in '",
      site_col, "':\n",
      paste(missing_sites, collapse = ", ")
    )
  }

  site_data_filtered <- if (!is.null(group_vars)) {
    site_data %>% dplyr::select(dplyr::all_of(c(site_col, group_vars)))
  } else {
    site_data
  }

  tibble::as_tibble(richness_matrix) %>%
    dplyr::mutate(!!site_col := rownames(richness_matrix)) %>%
    tidyr::pivot_longer(
      cols      = -dplyr::all_of(site_col),
      names_to  = "sample",
      values_to = "richness"
    ) %>%
    dplyr::left_join(site_data_filtered, by = site_col) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(site_col, group_vars)))) %>%
    dplyr::summarise(
      mean   = mean(richness),
      sd     = stats::sd(richness),
      median = stats::median(richness),
      cil    = stats::quantile(richness, 0.025),
      ciu    = stats::quantile(richness, 0.975),
      q05    = stats::quantile(richness, 0.05),
      q95    = stats::quantile(richness, 0.95),
      .groups = "drop"
    )
}


#' Build long-format richness data for density plotting
#'
#' Unlike \code{build_richness_summary()}, which summarises richness per site,
#' this function retains all posterior MCMC samples and joins site-level
#' metadata, producing a long-format tibble suitable for ridgeline or density
#' plots.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta") and Thiago Bicudo
#'
#' @param richness_matrix A \code{[sites x MCMC]} numeric matrix, output from
#'   \code{extract_richness()}.
#' @param site_data A data frame with site-level metadata. Must include a
#'   column matching \code{site_col} and all columns named in \code{group_vars}.
#' @param group_vars Character vector of grouping variable names.
#' @param site_col Character. Name of the site identifier column.
#'   Defaults to \code{"sites"}.
#' @param aggregate Logical. If \code{FALSE} (default), all site-level MCMC
#'   samples are retained — the density reflects both ecological variation
#'   among sites within a group and posterior uncertainty. If \code{TRUE},
#'   richness is averaged across sites within each group per MCMC sample,
#'   producing a density of the group mean only (within-group site variation
#'   is collapsed). \code{sum} is intentionally avoided because regions may
#'   differ in site count, which would confound totals with sampling effort.
#'
#' @return A long-format \code{tibble} with columns for the grouping variables,
#'   \code{sample} (MCMC iteration identifier), and \code{richness}.
#'
#' @details
#' Two modes are available via \code{aggregate}:
#'
#' \strong{aggregate = FALSE}: each row is one site x one MCMC sample. The
#' resulting density captures both ecological variation among sites within a
#' group and posterior uncertainty. Multiple modes in the density reflect real
#' within-group heterogeneity. Recommended for ridgeline plots.
#'
#' \strong{aggregate = TRUE}: each row is one group x one MCMC sample.
#' The density captures only posterior uncertainty around the group mean.
#' Useful when the question is "what is the expected mean richness of a site
#' in this region?" rather than the distribution of individual sites.
#'
#' @seealso \code{\link{extract_richness}}, \code{\link{plot_richness_ridgeline}}
#'
#' @export
build_richness_density <- function(richness_matrix,
                                   site_data,
                                   group_vars,
                                   site_col  = "sites",
                                   aggregate = FALSE) {

  stopifnot(
    "richness_matrix must be a numeric matrix"     = is.matrix(richness_matrix),
    "site_data must be a data frame"               = is.data.frame(site_data),
    "group_vars must be a character vector"        = is.character(group_vars) &&
      length(group_vars) >= 1,
    "site_col must be a single string"             = is.character(site_col) &&
      length(site_col) == 1,
    "site_col must exist in site_data"             = site_col %in% names(site_data),
    "aggregate must be logical"                    = is.logical(aggregate)
  )

  missing_sites <- setdiff(rownames(richness_matrix), site_data[[site_col]])
  if (length(missing_sites) > 0) {
    stop(
      "Rownames not found in '", site_col, "':\n",
      paste(missing_sites, collapse = ", ")
    )
  }

  missing_vars <- setdiff(group_vars, names(site_data))
  if (length(missing_vars) > 0) {
    stop(
      "The following group_vars columns are missing from site_data:\n",
      paste(missing_vars, collapse = ", ")
    )
  }

  richness_long <- tibble::as_tibble(richness_matrix) %>%
    dplyr::mutate(!!site_col := rownames(richness_matrix)) %>%
    tidyr::pivot_longer(
      cols      = -dplyr::all_of(site_col),
      names_to  = "sample",
      values_to = "richness"
    ) %>%
    dplyr::left_join(
      site_data %>% dplyr::select(dplyr::all_of(c(site_col, group_vars))),
      by = site_col
    )

  if (aggregate) {
    richness_long <- richness_long %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c("sample", group_vars)))) %>%
      dplyr::summarise(richness = mean(richness), .groups = "drop")
  }

  richness_long
}


#' Save a ggplot object with consistent settings
#'
#' Saves a \code{ggplot} object to a specified path using consistent output
#' dimensions, resolution, and background settings. Standardizes figure
#' saving across analyses.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta") and Thiago Bicudo
#'
#' @param plot_object A \code{ggplot} object to be saved.
#' @param filename Character. Name of the output file including extension
#'   (e.g., \code{"richness_plot.png"}).
#' @param output_dir Character. Directory where the file will be saved.
#'   Defaults to \code{getwd()}.
#' @param width Numeric. Plot width. Defaults to \code{20}.
#' @param height Numeric. Plot height. Defaults to \code{15}.
#' @param units Character. Units for width and height. One of \code{"cm"},
#'   \code{"in"}, or \code{"mm"}. Defaults to \code{"cm"}.
#' @param dpi Numeric. Resolution in dots per inch. Defaults to \code{600}.
#'
#' @return Invisibly returns the full file path to the saved plot.
#'
#' @seealso \code{\link[ggplot2]{ggsave}}
#'
#' @importFrom ggplot2 ggsave
#'
#' @export
save_plot <- function(plot_object,
                      filename,
                      output_dir = getwd(),
                      width      = 20,
                      height     = 15,
                      units      = "cm",
                      dpi        = 600) {

  stopifnot(
    "plot_object must be a ggplot"     = inherits(plot_object, "ggplot"),
    "filename must be a single string" = is.character(filename) &&
      length(filename) == 1,
    "output_dir must be a single string" = is.character(output_dir) &&
      length(output_dir) == 1,
    "dpi must be a positive number"    = is.numeric(dpi) && dpi > 0
  )

  if (!endsWith(output_dir, "/")) output_dir <- paste0(output_dir, "/")

  full_path <- normalizePath(
    paste0(output_dir, filename),
    winslash  = "/",
    mustWork  = FALSE
  )

  ggplot2::ggsave(
    filename = full_path,
    plot     = plot_object,
    width    = width,
    height   = height,
    units    = units,
    dpi      = dpi,
    bg       = "white"
  )

  message("Plot saved: ", full_path)
  invisible(full_path)
}
