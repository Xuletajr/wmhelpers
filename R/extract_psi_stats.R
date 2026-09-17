#' Extract Site-Level Occupancy Probability Statistics
#'
#' Computes posterior summaries (mean, SD, and 95\% credible interval) from
#' the \code{psi.samples} array of a fitted multi-species occupancy model.
#' Returns one matrix per statistic, with species as rows and sites as columns.
#'
#' @author Thiago Bicudo
#'
#' @param out_bestmodel Fitted model object containing \code{psi.samples}
#'   (a 3-dimensional array of dimensions
#'   \code{[samples x species x sites]}) and \code{sp.names}.
#' @param site_names Character vector of site names to use as column labels.
#'   Must match the number of sites in \code{psi.samples} (third dimension).
#'
#' @return A named list with four matrices, each of dimensions
#'   \code{[species x sites]} with species names as row names and
#'   \code{site_names} as column names:
#'   \code{mean} (posterior mean), \code{sd} (posterior standard deviation),
#'   \code{cri025} (lower bound of the 95\% credible interval), and
#'   \code{cri975} (upper bound of the 95\% credible interval).
#'
#' @details
#' This function operates on \code{psi.samples}, which represents site-level
#' occupancy probability and is shared across detection methods in integrated
#' models. It is not affected by the block structure of \code{alpha.samples}
#' and requires no method-specific slicing.
#'
#' Summary statistics are computed across the first dimension (MCMC samples)
#' of the 3D array using \code{apply()}.
#'
#' @importFrom stats sd quantile
#'
#' @export
# Extract psi stats
extract_psi_stats <- function(out_bestmodel, site_names, year_names = NULL) {
  stopifnot(
    `Model must contain psi.samples` = "psi.samples" %in% names(out_bestmodel),
    `Model must contain sp.names` = "sp.names" %in% names(out_bestmodel),
    `site_names must be character` = is.character(site_names)
  )

  psi_samples <- out_bestmodel$psi.samples
  n_dims <- length(dim(psi_samples))

  if (!n_dims %in% c(3, 4)) {
    stop("psi.samples must be 3D (samples x species x sites, single-season) ",
         "or 4D (samples x species x sites x years, multi-season).")
  }

  n_sites <- dim(psi_samples)[3]
  sp_names <- out_bestmodel$sp.names
  if (length(site_names) != n_sites) {
    stop("Length of site_names must match number of sites in psi.samples")
  }

  if (n_dims == 3) {
    # Single-season: original behaviour, byte-identical output.
    margin <- c(2, 3)
    dn <- list(sp_names, site_names)
  } else {
    # Multi-season: keep the year dimension in every summary statistic too.
    n_years <- dim(psi_samples)[4]
    if (is.null(year_names)) {
      stop("psi.samples is 4D (multi-season) -- `year_names` must be supplied ",
           "(e.g. dimnames(data.list$y)[[3]]).")
    }
    if (length(year_names) != n_years) {
      stop("Length of year_names must match number of years in psi.samples")
    }
    margin <- c(2, 3, 4)
    dn <- list(sp_names, site_names, year_names)
  }

  psi_mean   <- apply(psi_samples, margin, mean, na.rm = TRUE)
  psi_sd     <- apply(psi_samples, margin, sd, na.rm = TRUE)
  psi_cri025 <- apply(psi_samples, margin, quantile, probs = 0.025, na.rm = TRUE)
  psi_cri975 <- apply(psi_samples, margin, quantile, probs = 0.975, na.rm = TRUE)

  dimnames(psi_mean) <- dn
  dimnames(psi_sd) <- dn
  dimnames(psi_cri025) <- dn
  dimnames(psi_cri975) <- dn

  list(mean = psi_mean, sd = psi_sd, cri025 = psi_cri025, cri975 = psi_cri975)
}
