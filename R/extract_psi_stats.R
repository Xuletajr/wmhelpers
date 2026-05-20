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
extract_psi_stats <- function(out_bestmodel, site_names) {

  # ── Input validation ───────────────────────────────────────────────
  stopifnot(
    "Model must contain psi.samples" = "psi.samples" %in% names(out_bestmodel),
    "Model must contain sp.names"    = "sp.names" %in% names(out_bestmodel),
    "site_names must be character"   = is.character(site_names)
  )

  psi_samples <- out_bestmodel$psi.samples

  # Check dimensions
  if (length(dim(psi_samples)) != 3) {
    stop("psi.samples must be a 3D array: samples × species × sites")
  }

  n_species <- dim(psi_samples)[2]
  n_sites   <- dim(psi_samples)[3]

  if (length(site_names) != n_sites) {
    stop("Length of site_names must match number of sites in psi.samples")
  }

  # ── Compute statistics ─────────────────────────────────────────────
  psi_mean <- apply(psi_samples, c(2, 3), mean, na.rm = TRUE)
  psi_sd   <- apply(psi_samples, c(2, 3), sd,   na.rm = TRUE)

  psi_cri025 <- apply(psi_samples, c(2, 3), quantile, probs = 0.025, na.rm = TRUE)
  psi_cri975 <- apply(psi_samples, c(2, 3), quantile, probs = 0.975, na.rm = TRUE)

  # ── Assign dimnames ────────────────────────────────────────────────
  sp_names <- out_bestmodel$sp.names

  dimnames(psi_mean)   <- list(sp_names, site_names)
  dimnames(psi_sd)     <- list(sp_names, site_names)
  dimnames(psi_cri025) <- list(sp_names, site_names)
  dimnames(psi_cri975) <- list(sp_names, site_names)

  # ── Return as list ─────────────────────────────────────────────────
  return(list(
    mean   = psi_mean,
    sd     = psi_sd,
    cri025 = psi_cri025,
    cri975 = psi_cri975
  ))
}
