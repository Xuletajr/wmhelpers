#' Extract Site-Level Occupancy Probability Statistics
#'
#' Computes summary statistics (mean, SD, and 95\% credible intervals)
#' for site-level occupancy probabilities from posterior samples.
#' Extracts posterior samples from the 3D array
#' (\code{samples × species × sites}) and computes summary statistics
#' across the MCMC dimension.
#'
#' @param out_bestmodel Fitted model object containing
#'   \code{psi.samples} and \code{sp.names}.
#' @param site_names Character vector of site names used as column labels.
#'
#' @return Named list with four matrices:
#' \describe{
#'   \item{\code{mean}}{Posterior mean occupancy probabilities}
#'   \item{\code{sd}}{Posterior standard deviation}
#'   \item{\code{cri025}}{Lower bound of 95\% credible interval}
#'   \item{\code{cri975}}{Upper bound of 95\% credible interval}
#' }
#' Each matrix has species as rows and sites as columns.
#'
#' @author Thiago Bicudo
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
