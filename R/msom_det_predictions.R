#' Predict Community-Level Detection Probability Across a Continuous
#' Covariate and an Optional Categorical Class
#'
#' Analogous to \code{\link{compute_occupancy}}, but for the detection
#' sub-model (\code{alpha.comm.samples} / \code{X.p}) instead of occupancy
#' (\code{beta.comm.samples} / \code{X}). Supports both single-method models
#' (\code{msPGOcc}, where \code{X.p} is a single matrix) and integrated
#' multi-method models (\code{intMsPGOcc}, where \code{X.p} is a list of
#' matrices, one per method).
#'
#' @author Thiago Bicudo
#'
#' @param out_bestmodel Fitted model object containing
#'   \code{alpha.comm.samples} and \code{X.p}.
#' @param scaled_centered_list Named list of scaled covariate sequences, as
#'   returned by \code{\link{build_occ_sequences}}.
#' @param variables_to_scale Character vector of continuous variable names.
#'   Each name must appear in both \code{scaled_centered_list} and the column
#'   names of the relevant \code{alpha.comm.samples} block.
#' @param class_levels_labels Character vector of display labels for the
#'   categorical detection covariate (e.g. \code{c("Dry season", "Wet
#'   season")}). The first element is treated as the reference/intercept
#'   level. Must be parallel to \code{class_levels_model}.
#' @param class_levels_model Character vector of coefficient names exactly as
#'   they appear in \code{alpha.comm.samples} for the target method block
#'   (e.g. \code{c("(Intercept)", "season")}). Must be parallel to
#'   \code{class_levels_labels}.
#' @param nsamp Integer. Number of posterior samples to use for prediction.
#' @param n_pred Integer. Number of prediction points per variable.
#'   Defaults to \code{100L}.
#' @param method_index Integer. Index of the detection method block in
#'   \code{alpha.comm.samples}/\code{X.p} (1 = first method). Only relevant
#'   for integrated multi-method models; must be \code{1L} for single-method
#'   models. Defaults to \code{1L}.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\code{predDet}}{3-dimensional array \code{[n_pred x nsamp x
#'     n_combinations]} of posterior detection probabilities on the
#'     probability scale (after \code{plogis}).}
#'   \item{\code{comb_df}}{Data frame of all variable x class_label
#'     combinations used to build the array.}
#' }
#'
#' @seealso \code{\link{build_occ_sequences}}, \code{\link{compute_occupancy}},
#'   \code{\link{build_occ_summary}}, \code{\link{compute_det_predictions}},
#'   \code{\link{predict_detection_julian}}
#'
#' @importFrom stats plogis
#'
#' @export
compute_detection <- function(out_bestmodel,
                              scaled_centered_list,
                              variables_to_scale,
                              class_levels_labels,
                              class_levels_model,
                              nsamp,
                              n_pred       = 100L,
                              method_index = 1L) {

  # ── Input validation ─────────────────────────────────────────────────────────
  stopifnot(
    "out_bestmodel must contain alpha.comm.samples" =
      !is.null(out_bestmodel$alpha.comm.samples),
    "out_bestmodel must contain X.p"                =
      !is.null(out_bestmodel$X.p),
    "scaled_centered_list must be a list"           = is.list(scaled_centered_list),
    "variables_to_scale must be character"          = is.character(variables_to_scale),
    "class_levels_labels must be character"         = is.character(class_levels_labels),
    "class_levels_model must be character"          = is.character(class_levels_model),
    "class_levels_labels and class_levels_model must have the same length" =
      length(class_levels_labels) == length(class_levels_model),
    "nsamp must be numeric"                         = is.numeric(nsamp),
    "n_pred must be numeric"                        = is.numeric(n_pred)
  )

  # ── Select community samples for the target method ──────────────────────────
  # Integrated models (intMsPGOcc): X.p is a list of matrices, one per method,
  # and alpha.comm.samples has one column block per method (with possibly
  # duplicated names, e.g. two "(Intercept)"), so blocks must be selected by
  # cumulative position.
  # Single-method models (msPGOcc): X.p is a single matrix, and
  # alpha.comm.samples already contains only that method's coefficients —
  # there are no blocks to select.
  alpha_comm <- out_bestmodel$alpha.comm.samples

  if (is.list(out_bestmodel$X.p) && !is.data.frame(out_bestmodel$X.p)) {

    block_sizes  <- sapply(out_bestmodel$X.p, ncol)
    block_starts <- c(1, cumsum(block_sizes[-length(block_sizes)]) + 1)

    if (method_index < 1L || method_index > length(out_bestmodel$X.p)) {
      stop(
        "method_index must be a valid method position (1 to ",
        length(out_bestmodel$X.p), ")."
      )
    }

    col_start <- block_starts[method_index]
    col_end   <- col_start + block_sizes[method_index] - 1

    alpha_samples <- alpha_comm[, col_start:col_end, drop = FALSE]

  } else {

    if (method_index != 1L) {
      stop(
        "method_index > 1 is only valid for integrated multi-method models ",
        "(where out_bestmodel$X.p is a list). This model has a single X.p matrix."
      )
    }

    alpha_samples <- alpha_comm
  }

  # ── Validate nsamp ───────────────────────────────────────────────────────────
  if (nsamp > nrow(alpha_samples)) {
    stop(
      "nsamp (", nsamp, ") cannot exceed number of posterior samples (",
      nrow(alpha_samples), ")."
    )
  }

  alpha_samples <- alpha_samples[seq_len(nsamp), , drop = FALSE]
  coef_names    <- colnames(alpha_samples)

  # ── Build combinations: variable x class_label ───────────────────────────────
  comb_df <- expand.grid(
    variable   = variables_to_scale,
    year_label = class_levels_labels,   # generic slot, feeds build_occ_summary()
    stringsAsFactors = FALSE
  )

  class_label_to_model <- stats::setNames(class_levels_model, class_levels_labels)
  n_comb               <- nrow(comb_df)

  predDet <- array(
    NA_real_,
    dim      = c(n_pred, nsamp, n_comb),
    dimnames = list(
      NULL,
      NULL,
      apply(comb_df, 1, function(x) paste(x, collapse = "_"))
    )
  )

  # ── Prediction loop ───────────────────────────────────────────────────────────
  for (c in seq_len(n_comb)) {

    var         <- comb_df$variable[c]
    class_label <- comb_df$year_label[c]
    class_coef  <- class_label_to_model[class_label]
    x           <- scaled_centered_list[[var]]

    for (i in seq_len(nsamp)) {

      eta <- alpha_samples[i, "(Intercept)"]

      # Continuous variable effect
      if (var %in% coef_names) {
        eta <- eta + alpha_samples[i, var] * x
      }

      # Class effect — skip if reference level
      if (!is.null(class_coef) &&
          class_coef != "(Intercept)" &&
          class_coef %in% coef_names) {
        eta <- eta + alpha_samples[i, class_coef]
      }

      predDet[, i, c] <- stats::plogis(eta)
    }
  }

  list(predDet = predDet, comb_df = comb_df)
}
