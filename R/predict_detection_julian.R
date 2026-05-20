#' Predict Community-Level Detection Probability Across Julian Days
#'
#' Generates predicted detection probabilities across a sequence of Julian days
#' using community-level detection coefficients (linear and optional quadratic)
#' from a fitted occupancy model. Supports both single-method and integrated
#' multi-method models.
#'
#' @author Thiago Bicudo
#'
#' @param out_bestmodel Fitted model object containing
#'   \code{alpha.comm.samples} and \code{X.p}.
#' @param julian_day Matrix of Julian day values (sites x surveys). Used only
#'   to define the prediction range (min-max).
#' @param m_julian Numeric. Mean Julian day used for centering.
#' @param sd_julian Numeric. Standard deviation of Julian days used for scaling.
#' @param linear_coef Character. Name of the linear coefficient within the
#'   target method block. Default is \code{"days"}.
#' @param quadratic_coef Character or \code{NULL}. Name of the quadratic term
#'   within the target method block. Use \code{NULL} to omit.
#'   Default is \code{"I(days^2)"}.
#' @param n_pred Integer. Number of prediction points. Default is \code{100L}.
#' @param nsamp Integer. Number of posterior samples to use.
#' @param method_index Integer. Index of the detection method in the integrated
#'   model (1 = first method, 2 = second, etc.). Used to select the correct
#'   column block from \code{alpha.comm.samples} by position, avoiding
#'   ambiguity caused by duplicate column names across methods (e.g., two
#'   \code{"(Intercept)"} columns). For single-method models, the default
#'   \code{1L} requires no change.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{pred_det_df}}{Data frame with columns:
#'     \code{pred_julian} (scaled), \code{seq_julian} (raw),
#'     \code{mean}, \code{criL}, \code{criU}.}
#'   \item{\code{seq_julian}}{Sequence of raw Julian days.}
#'   \item{\code{pred_julian}}{Scaled Julian day values.}
#' }
#'
#' @details
#' Detection probability is computed using the inverse-logit transformation
#' applied to the linear predictor:
#'
#' \deqn{logit(p) = \alpha_0 + \alpha_1 x + \alpha_2 x^2}
#'
#' where \eqn{x} is the standardized Julian day.
#'
#' In integrated models, \code{alpha.comm.samples} contains one column block
#' per detection method. Column names are frequently duplicated across methods
#' (e.g., two \code{"(Intercept)"} columns — one per method). Indexing by name
#' in such cases would silently return the wrong columns. This function selects
#' the correct block by cumulative position using \code{X.p[[method_index]]}
#' to determine each block's size, making name-based lookup unnecessary and
#' safe.
#'
#' For integrated models, the recommended calls are:
#' \preformatted{
#' # Audio (method_index = 1)
#' det_results_aud <- predict_detection_julian(
#'   out_bestmodel  = out_bestmodel,
#'   julian_day     = day_index_aud[, -1],
#'   m_julian       = m_day_aud,
#'   sd_julian      = sd_day_aud,
#'   linear_coef    = "days",
#'   quadratic_coef = "I(days^2)",
#'   n_pred         = 100,
#'   nsamp          = nsamp,
#'   method_index   = 1L
#' )
#'
#' # Camera trap (method_index = 2)
#' det_results_cam <- predict_detection_julian(
#'   out_bestmodel  = out_bestmodel,
#'   julian_day     = day_index_cam[, -1],
#'   m_julian       = m_day_cam,
#'   sd_julian      = sd_day_cam,
#'   linear_coef    = "numD",
#'   quadratic_coef = NULL,
#'   n_pred         = 100,
#'   nsamp          = nsamp,
#'   method_index   = 2L
#' )
#' }
#'
#' @importFrom stats plogis quantile
#'
#' @export
predict_detection_julian <- function(
    out_bestmodel,
    julian_day,
    m_julian,
    sd_julian,
    linear_coef    = "days",
    quadratic_coef = "I(days^2)",
    n_pred         = 100L,
    nsamp,
    method_index   = 1L
) {

  # ── Input validation ───────────────────────────────────────────────────────
  stopifnot(
    "out_bestmodel must contain alpha.comm.samples" =
      !is.null(out_bestmodel$alpha.comm.samples),
    "out_bestmodel must contain X.p"                =
      !is.null(out_bestmodel$X.p),
    "julian_day must be a matrix"                   = is.matrix(julian_day),
    "m_julian must be numeric"                      = is.numeric(m_julian),
    "sd_julian must be numeric"                     = is.numeric(sd_julian),
    "n_pred must be numeric"                        = is.numeric(n_pred),
    "nsamp must be numeric"                         = is.numeric(nsamp),
    "method_index must be a single integer"         = is.numeric(method_index) &&
      length(method_index) == 1,
    "method_index must be a valid method position"  =
      method_index >= 1L && method_index <= length(out_bestmodel$X.p)
  )

  # ── Select community samples for the target method ─────────────────────────
  # alpha.comm.samples contains one column block per detection method.
  # Column names are duplicated across methods (e.g., two "(Intercept)").
  # Block boundaries are computed from cumulative block sizes via X.p to
  # avoid silent errors from name-based indexing on duplicate column names.
  alpha_comm <- out_bestmodel$alpha.comm.samples

  block_sizes  <- sapply(out_bestmodel$X.p, ncol)
  block_starts <- c(1, cumsum(block_sizes[-length(block_sizes)]) + 1)

  col_start <- block_starts[method_index]
  col_end   <- col_start + block_sizes[method_index] - 1

  alpha_samples <- alpha_comm[, col_start:col_end, drop = FALSE]

  # ── Validate nsamp ─────────────────────────────────────────────────────────
  if (nsamp > nrow(alpha_samples)) {
    stop(
      "nsamp (", nsamp, ") cannot be greater than the number of ",
      "posterior samples (", nrow(alpha_samples), ")."
    )
  }

  # ── Validate coefficient names within the selected block ───────────────────
  coef_names <- colnames(alpha_samples)

  if (!(linear_coef %in% coef_names)) {
    stop(
      "Linear coefficient '", linear_coef, "' not found in alpha.comm.samples ",
      "for method_index = ", method_index, ".\n",
      "Available in this block: ", paste(coef_names, collapse = ", ")
    )
  }

  if (!is.null(quadratic_coef) && !(quadratic_coef %in% coef_names)) {
    stop(
      "Quadratic coefficient '", quadratic_coef, "' not found in ",
      "alpha.comm.samples for method_index = ", method_index, ".\n",
      "Available in this block: ", paste(coef_names, collapse = ", ")
    )
  }

  # ── Subset posterior samples ───────────────────────────────────────────────
  alpha_samples <- alpha_samples[1:nsamp, , drop = FALSE]

  intercept <- alpha_samples[, "(Intercept)"]
  beta_lin  <- alpha_samples[, linear_coef]

  if (!is.null(quadratic_coef)) {
    beta_quad <- alpha_samples[, quadratic_coef]
  }

  # ── Create Julian day sequence ─────────────────────────────────────────────
  seq_julian <- seq(
    min(julian_day, na.rm = TRUE),
    max(julian_day, na.rm = TRUE),
    length.out = n_pred
  )

  pred_julian <- (seq_julian - m_julian) / sd_julian

  # ── Predict detection probability ──────────────────────────────────────────
  pred_mat <- matrix(NA, nrow = n_pred, ncol = nsamp)

  for (i in seq_len(nsamp)) {
    eta <- intercept[i] + beta_lin[i] * pred_julian
    if (!is.null(quadratic_coef)) {
      eta <- eta + beta_quad[i] * (pred_julian^2)
    }
    pred_mat[, i] <- plogis(eta)
  }

  # ── Summarise posterior predictions ────────────────────────────────────────
  pred_det_df <- data.frame(
    pred_julian = pred_julian,
    seq_julian  = seq_julian,
    mean        = apply(pred_mat, 1, mean),
    criL        = apply(pred_mat, 1, quantile, probs = 0.025),
    criU        = apply(pred_mat, 1, quantile, probs = 0.975)
  )

  return(list(
    pred_det_df = pred_det_df,
    seq_julian  = seq_julian,
    pred_julian = pred_julian
  ))
}
