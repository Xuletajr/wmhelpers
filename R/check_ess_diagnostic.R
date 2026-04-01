#' Check Effective Sample Size (ESS)
#'
#' Identifies parameters with ESS below a specified threshold and optionally
#' exports a CSV of problematic parameters.
#'
#' @author Thiago Bicudo
#'
#' @param ess_values Named numeric vector of ESS values
#'   (e.g. \code{out_bestmodel$ESS$beta}).
#'
#' @param parameter_name Character label for messages (e.g. \code{"BETA"}).
#'
#' @param threshold Numeric. ESS values below this are flagged.
#'   Defaults to \code{450} (150 per chain × 3 chains).
#'
#' @param output_path Character. If provided, the problematic table is written
#'   to this CSV path. Defaults to \code{NULL}.
#'
#' @return A \code{data.frame} with columns \code{Parameter} and \code{ESS}
#'   for all flagged parameters, or \code{NULL} if none found.
#'
#' @export
check_ess <- function(
    ess_values,
    parameter_name,
    threshold = 450,
    output_path = NULL
) {

  # ── Input validation ────────────────────────────────────────────────────────
  stopifnot(
    is.numeric(ess_values),
    !is.null(names(ess_values)),
    is.character(parameter_name),
    length(parameter_name) == 1
  )

  # ── Identify problematic parameters ─────────────────────────────────────────
  low_ess <- ess_values[ess_values < threshold]

  if (length(low_ess) == 0) {
    message(paste0("All ", parameter_name, " parameters have ESS ≥ ", threshold))
    return(NULL)
  }

  # ── Build output table ──────────────────────────────────────────────────────
  result <- data.frame(
    Parameter = names(low_ess),
    ESS = as.numeric(low_ess),
    row.names = NULL
  )

  # ── Warning message ─────────────────────────────────────────────────────────
  warning(
    paste0(
      "Some ", parameter_name,
      " parameters have ESS < ", threshold,
      ". Consider increasing iterations or improving mixing."
    )
  )

  # ── Optional export ─────────────────────────────────────────────────────────
  if (!is.null(output_path)) {
    utils::write.csv(result, output_path, row.names = FALSE)
  }

  return(result)
}
