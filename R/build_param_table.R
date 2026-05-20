#' Build Species-Level Parameter Summary Table
#'
#' Creates a tidy summary table (mean, SD, 95\% credible interval, and effect
#' direction) for species-level MCMC samples from occupancy (beta) or detection
#' (alpha) parameters.
#'
#' For integrated multi-method models, pass the output of
#' \code{extract_method_samples()} rather than the full \code{alpha.samples},
#' so that \code{n_species} and \code{n_covariates} correspond to a single
#' method with a uniform block structure.
#'
#' @author Thiago Bicudo
#'
#' @param samples A matrix or \code{mcmc} object of MCMC samples (e.g.
#'   \code{out_bestmodel$beta.samples} or a slice of
#'   \code{out_bestmodel$alpha.samples} from \code{extract_method_samples()}).
#' @param sp_names Character vector of species names corresponding to the
#'   method being summarised. For occupancy, use \code{out_bestmodel$sp.names}.
#'   For detection in integrated models, use the species vector for the
#'   relevant method (e.g. \code{rownames(data.list$y$audio)}).
#' @param n_covariates Integer. Total number of parameters per species in the
#'   design matrix, including the intercept
#'   (\code{ncol(out_bestmodel$X)} for beta;
#'   \code{ncol(out_bestmodel$X.p[[i]])} for alpha of method \code{i}).
#' @param output_path Character. If provided, the table is exported as a CSV.
#'   Defaults to \code{NULL}.
#'
#' @return A \code{data.frame} with columns:
#'   \code{param_spp_code}, \code{species}, \code{Mean}, \code{SD},
#'   \code{CRI_low}, \code{CRI_high}, \code{Effect}.
#'
#' @details
#' The first \code{length(sp_names)} columns of \code{samples} are assumed to
#' correspond to intercepts, for which \code{Effect} is set to \code{""}.
#' For all subsequent parameter blocks, \code{Effect} is classified as
#' \code{"Positive"}, \code{"Negative"}, or \code{"No effect"} depending on
#' whether the 95\% credible interval excludes zero.
#'
#' The recommended workflow for integrated models is:
#' \preformatted{
#' # Audio
#' build_param_table(
#'   samples      = extract_method_samples(alpha.samples, method_spec, "audio"),
#'   sp_names     = rownames(data.list$y$audio),
#'   n_covariates = ncol(out_bestmodel$X.p[[1]])
#' )
#'
#' # Camera trap
#' build_param_table(
#'   samples      = extract_method_samples(alpha.samples, method_spec, "camtrap"),
#'   sp_names     = rownames(data.list$y$camtrap),
#'   n_covariates = ncol(out_bestmodel$X.p[[2]])
#' )
#' }
#'
#' @importFrom stats quantile sd
#'
#' @export
build_param_table <- function(samples,
                               sp_names,
                               n_covariates,
                               output_path = NULL) {

  # ── Input validation ───────────────────────────────────────────────────────
  if (inherits(samples, "mcmc")) samples <- as.matrix(samples)

  stopifnot(
    "samples must be a matrix or mcmc object"  = is.matrix(samples),
    "sp_names must be a character vector"      = is.character(sp_names),
    "n_covariates must be a positive integer"  = is.numeric(n_covariates) &&
      n_covariates >= 1
  )

  n_species <- length(sp_names)
  n_params  <- n_species * n_covariates

  if (ncol(samples) != n_params) {
    stop(
      "ncol(samples) = ", ncol(samples), " does not match ",
      "length(sp_names) x n_covariates = ", n_species, " x ", n_covariates,
      " = ", n_params, ".\n",
      "For integrated models, pass a method slice from extract_method_samples()."
    )
  }

  # ── Summary statistics ─────────────────────────────────────────────────────
  mean_vals <- apply(samples, 2, mean)
  sd_vals   <- apply(samples, 2, sd)
  cri_low   <- apply(samples, 2, quantile, probs = 0.025)
  cri_high  <- apply(samples, 2, quantile, probs = 0.975)

  # ── Column names and species vector ───────────────────────────────────────
  param_spp_code <- colnames(samples)

  if (is.null(param_spp_code)) {
    param_spp_code <- paste0("param_", seq_len(ncol(samples)))
  }

  species_rep <- rep(sp_names, times = n_covariates)

  # ── Effect direction ───────────────────────────────────────────────────────
  effect <- ifelse(
    cri_low  > 0, "Positive",
    ifelse(cri_high < 0, "Negative", "No effect")
  )

  # First block (intercepts) gets no effect label
  effect[seq_len(n_species)] <- ""

  # ── Assemble output ────────────────────────────────────────────────────────
  out <- data.frame(
    param_spp_code   = param_spp_code,
    species          = species_rep,
    Mean             = round(mean_vals, 3),
    SD               = round(sd_vals,   3),
    CRI_low          = round(cri_low,   3),
    CRI_high         = round(cri_high,  3),
    Effect           = effect,
    stringsAsFactors = FALSE
  )

  # ── Export ─────────────────────────────────────────────────────────────────
  if (!is.null(output_path)) {
    utils::write.csv(out, output_path, row.names = FALSE)
    message("Parameter table saved: ", basename(output_path))
  }

  return(out)
}
