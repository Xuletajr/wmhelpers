#' Predict Community-Level Detection Probability Across Julian Days
#'
#' Generates predicted detection probabilities across a sequence of Julian days
#' using community-level detection coefficients (linear and optional quadratic)
#' from a fitted occupancy model.
#'
#' @author Thiago Bicudo
#'
#' @param out_bestmodel Fitted model object containing
#'   \code{alpha.comm.samples}.
#'
#' @param julian_day Matrix of Julian day values (sites × surveys). Used only
#'   to define the prediction range (min–max).
#'
#' @param m_julian Numeric. Mean Julian day used for centering.
#'
#' @param sd_julian Numeric. Standard deviation of Julian days used for scaling.
#'
#' @param linear_coef Character. Name of the linear Julian day coefficient in
#'   \code{alpha.comm.samples}. Default is \code{"julian_std"}.
#'
#' @param quadratic_coef Character or \code{NULL}. Name of the quadratic term.
#'   Use \code{NULL} to omit. Default is \code{"I(julian_std^2)"}.
#'
#' @param n_pred Integer. Number of prediction points. Default is \code{100L}.
#'
#' @param nsamp Integer. Number of posterior samples to use.
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
#' \deqn{
#' logit(p) = \alpha_0 + \alpha_1 x + \alpha_2 x^2
#' }
#'
#' where \eqn{x} is the standardized Julian day.
#'
#' @importFrom stats plogis quantile
#'
#' @export
predict_detection_julian <- function(
    out_bestmodel,
    julian_day,
    m_julian,
    sd_julian,
    linear_coef = "julian_std",
    quadratic_coef = "I(julian_std^2)",
    n_pred = 100L,
    nsamp
) {

  # ── Input validation ────────────────────────────────────────────────────────
  stopifnot(
    !is.null(out_bestmodel$alpha.comm.samples),
    is.matrix(julian_day),
    is.numeric(m_julian),
    is.numeric(sd_julian),
    is.numeric(n_pred),
    is.numeric(nsamp)
  )

  alpha_samples <- out_bestmodel$alpha.comm.samples

  if (nsamp > nrow(alpha_samples)) {
    stop("nsamp cannot be greater than number of posterior samples")
  }

  # ── Extract coefficients ────────────────────────────────────────────────────
  coef_names <- colnames(alpha_samples)

  if (!(linear_coef %in% coef_names)) {
    stop("Linear coefficient not found in alpha.comm.samples")
  }

  if (!is.null(quadratic_coef) && !(quadratic_coef %in% coef_names)) {
    stop("Quadratic coefficient not found in alpha.comm.samples")
  }

  # Subset samples
  alpha_samples <- alpha_samples[1:nsamp, , drop = FALSE]

  intercept <- alpha_samples[, "(Intercept)"]
  beta_lin  <- alpha_samples[, linear_coef]

  if (!is.null(quadratic_coef)) {
    beta_quad <- alpha_samples[, quadratic_coef]
  }

  # ── Create Julian sequence ──────────────────────────────────────────────────
  seq_julian <- seq(
    min(julian_day, na.rm = TRUE),
    max(julian_day, na.rm = TRUE),
    length.out = n_pred
  )

  # Standardize
  pred_julian <- (seq_julian - m_julian) / sd_julian

  # ── Predict detection probability ───────────────────────────────────────────
  pred_mat <- matrix(NA, nrow = n_pred, ncol = nsamp)

  for (i in seq_len(nsamp)) {

    eta <- intercept[i] + beta_lin[i] * pred_julian

    if (!is.null(quadratic_coef)) {
      eta <- eta + beta_quad[i] * (pred_julian^2)
    }

    pred_mat[, i] <- plogis(eta)
  }

  # ── Summarise posterior predictions ─────────────────────────────────────────
  mean_pred <- apply(pred_mat, 1, mean)
  criL      <- apply(pred_mat, 1, quantile, probs = 0.025)
  criU      <- apply(pred_mat, 1, quantile, probs = 0.975)

  pred_det_df <- data.frame(
    pred_julian = pred_julian,
    seq_julian  = seq_julian,
    mean        = mean_pred,
    criL        = criL,
    criU        = criU
  )

  # ── Return ──────────────────────────────────────────────────────────────────
  return(list(
    pred_det_df = pred_det_df,
    seq_julian  = seq_julian,
    pred_julian = pred_julian
  ))
}
