#' Check R-hat Convergence Diagnostic
#'
#' Evaluates R-hat (Gelman-Rubin) values for a set of MCMC parameters and
#' returns a table of any parameters exceeding the acceptable threshold.
#' Issues a \code{message()} when all values are acceptable and a
#' \code{warning()} when any exceed the threshold.
#'
#' @author Thiago Bicudo
#'
#' @param rhat_values Named numeric vector of R-hat values
#'   (e.g. \code{out_bestmodel$rhat$beta}).
#'
#' @param samples Matrix of MCMC samples whose column names identify
#'   parameters (e.g. \code{out_bestmodel$beta.samples}). Used only to
#'   recover parameter names.
#'
#' @param parameter_name Character label used in messages
#'   (e.g. \code{"BETA"} or \code{"ALPHA"}).
#'
#' @param threshold Numeric. R-hat values at or above this value are flagged.
#'   Defaults to \code{1.1}.
#'
#' @return A \code{data.frame} with columns \code{Parameter} and \code{Rhat}
#'   for all problematic parameters, or \code{NULL} if none are found.
#'
#' @export
check_rhat <- function(
    rhat_values,
    samples,
    parameter_name,
    threshold = 1.1
) {

  # ── Input validation ────────────────────────────────────────────────────────
  stopifnot(
    is.numeric(rhat_values),
    is.matrix(samples),
    is.character(parameter_name),
    length(parameter_name) == 1
  )

  # ── Recover parameter names ─────────────────────────────────────────────────
  param_names <- colnames(samples)

  if (is.null(param_names)) {
    stop("samples must have column names corresponding to parameters")
  }

  if (length(rhat_values) != length(param_names)) {
    stop("Length of rhat_values must match number of parameters (columns in samples)")
  }

  # ── Identify problematic parameters ─────────────────────────────────────────
  problematic_idx <- which(rhat_values >= threshold)

  # ── All good case ───────────────────────────────────────────────────────────
  if (length(problematic_idx) == 0) {
    message(
      sprintf("All %s parameters converged (R-hat < %.2f).",
              parameter_name, threshold)
    )
    return(NULL)
  }

  # ── Build output table ──────────────────────────────────────────────────────
  out <- data.frame(
    Parameter = param_names[problematic_idx],
    Rhat = round(rhat_values[problematic_idx], 3),
    row.names = NULL
  )

  # ── Warning ─────────────────────────────────────────────────────────────────
  warning(
    sprintf(
      "%d %s parameter(s) with R-hat >= %.2f detected.",
      nrow(out), parameter_name, threshold
    )
  )

  return(out)
}
