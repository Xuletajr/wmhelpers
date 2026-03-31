#' Build Species-Level Parameter Summary Table
#'
#' Creates a tidy summary table (mean, SD, 95\% credible interval, and effect
#' direction) for species-level MCMC samples from occupancy (beta) or detection
#' (alpha) parameters.
#'
#' @author Thiago Bicudo
#'
#' @param samples Matrix of MCMC samples (e.g.
#'   \code{out_bestmodel$beta.samples} or
#'   \code{out_bestmodel$alpha.samples}).
#'
#' @param sp_names Character vector of species names in the same order as
#'   \code{out_bestmodel$sp.names}.
#'
#' @param n_covariates Integer. Number of covariates in the design matrix
#'   (\code{ncol(out_bestmodel$X)} for beta;
#'   \code{ncol(out_bestmodel$X.p)} for alpha).
#'
#' @param output_path Character. If provided, the table is exported as a CSV.
#'   Defaults to \code{NULL}.
#'
#' @return A \code{data.frame} with columns:
#' \code{param_spp_code}, \code{species}, \code{Mean}, \code{SD},
#' \code{CRI_low}, \code{CRI_high}, \code{Effect}.
#'
#' @details
#' The first \code{length(sp_names)} rows correspond to intercepts, for which
#' \code{Effect} is set to \code{""}. For covariates, \code{Effect} is classified
#' as \code{"Positive"}, \code{"Negative"}, or \code{"No effect"} depending on
#' whether the 95\% credible interval excludes zero.
#'
#' @importFrom stats quantile sd
#'
#' @export
build_param_table <- function(
    samples,
    sp_names,
    n_covariates,
    output_path = NULL
) {

  # ── Input validation ────────────────────────────────────────────────────────
  stopifnot(
    is.matrix(samples),
    is.character(sp_names),
    is.numeric(n_covariates),
    n_covariates >= 1
  )

  n_species <- length(sp_names)
  n_params  <- n_species * n_covariates

  if (ncol(samples) != n_params) {
    stop("Number of columns in samples does not match species × covariates")
  }

  # ── Summary statistics ──────────────────────────────────────────────────────
  mean_vals <- apply(samples, 2, mean)
  sd_vals   <- apply(samples, 2, sd)
  cri_low   <- apply(samples, 2, quantile, probs = 0.025)
  cri_high  <- apply(samples, 2, quantile, probs = 0.975)

  # ── Build base table ────────────────────────────────────────────────────────
  param_spp_code <- colnames(samples)

  if (is.null(param_spp_code)) {
    param_spp_code <- paste0("param_", seq_len(ncol(samples)))
  }

  species_rep <- rep(sp_names, times = n_covariates)

  # ── Determine effect direction ──────────────────────────────────────────────
  effect <- ifelse(
    cri_low > 0, "Positive",
    ifelse(cri_high < 0, "Negative", "No effect")
  )

  # Intercepts → no effect label
  effect[1:n_species] <- ""

  # ── Final table ─────────────────────────────────────────────────────────────
  out <- data.frame(
    param_spp_code = param_spp_code,
    species        = species_rep,
    Mean           = mean_vals,
    SD             = sd_vals,
    CRI_low        = cri_low,
    CRI_high       = cri_high,
    Effect         = effect,
    stringsAsFactors = FALSE
  )

  # ── Export if requested ─────────────────────────────────────────────────────
  if (!is.null(output_path)) {
    utils::write.csv(out, output_path, row.names = FALSE)
  }

  return(out)
}
