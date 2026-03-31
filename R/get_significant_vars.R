#' Identify Significant Variables from Community-Level Beta Samples
#'
#' Returns variable names whose credible interval (CRI) excludes zero,
#' based on community-level beta posterior samples from occupancy models.
#' This is useful for selecting predictors for visualization and interpretation.
#'
#' @param beta_comm_samples Matrix of community-level beta MCMC samples
#'   (e.g., \code{out_bestmodel$beta.comm.samples}).
#'
#' @param probs Numeric vector of length 2 specifying the lower and upper
#'   quantile bounds for the credible interval. Defaults to \code{c(0.025, 0.975)}.
#'
#' @param exclude Character vector of column names to always exclude
#'   (e.g., intercept). Defaults to \code{"(Intercept)"}.
#'
#' @return Character vector of variable names with credible intervals
#'   that do not include zero.
#'
#' @details
#' A variable is considered "significant" if its posterior credible interval
#' (e.g., 95%) lies entirely above or below zero. This function operates at the
#' community level and does not account for species-specific variation.
#'
#' @examples
#' \dontrun{
#' sig_vars <- get_significant_vars(out_bestmodel$beta.comm.samples)
#' }
#'
#' @author Thiago Bicudo
#'
#' @export
get_significant_vars <- function(
    beta_comm_samples,
    probs = c(0.025, 0.975),
    exclude = "(Intercept)"
) {

  # Check input
  if (!is.matrix(beta_comm_samples)) {
    stop("beta_comm_samples must be a matrix.")
  }

  # Compute credible intervals
  ci <- apply(beta_comm_samples, 2, quantile, probs = probs, na.rm = TRUE)

  # Identify variables whose CI excludes zero
  significant <- (ci[1, ] > 0) | (ci[2, ] < 0)

  vars <- colnames(beta_comm_samples)[significant]

  # Remove excluded variables
  vars <- setdiff(vars, exclude)

  return(vars)
}
