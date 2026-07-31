# R/msom_det_predictions.R
# =============================================================================
# Detection prediction pipeline functions (mirrors msom_occ_predictions.R)
# =============================================================================

#' Predict Community-Level Detection Probability Across a Continuous
#' Covariate and an Optional Categorical Class
#'
#' Analogous to \code{\link{compute_occupancy}}, but for the detection
#' sub-model (\code{alpha.comm.samples} / \code{X.p}) instead of occupancy
#' (\code{beta.comm.samples} / \code{X}).
#'
#' @author Thiago Bicudo
#' @seealso \code{\link{build_occ_sequences}}, \code{\link{compute_occupancy}},
#'   \code{\link{build_occ_summary}}, \code{\link{compute_det_predictions}}
#' @importFrom stats plogis
#' @export
compute_detection <- function(out_bestmodel,
                              scaled_centered_list,
                              variables_to_scale,
                              class_levels_labels,
                              class_levels_model,
                              nsamp,
                              n_pred       = 100L,
                              method_index = 1L) {

  stopifnot(
    "out_bestmodel must contain alpha.comm.samples" = !is.null(out_bestmodel$alpha.comm.samples),
    "out_bestmodel must contain X.p"                = !is.null(out_bestmodel$X.p),
    "scaled_centered_list must be a list"            = is.list(scaled_centered_list),
    "variables_to_scale must be character"           = is.character(variables_to_scale),
    "class_levels_labels must be character"          = is.character(class_levels_labels),
    "class_levels_model must be character"           = is.character(class_levels_model),
    "class_levels_labels and class_levels_model must have the same length" =
      length(class_levels_labels) == length(class_levels_model),
    "nsamp must be numeric"                          = is.numeric(nsamp),
    "n_pred must be numeric"                         = is.numeric(n_pred)
  )

  alpha_comm   <- out_bestmodel$alpha.comm.samples
  block_sizes  <- sapply(out_bestmodel$X.p, ncol)
  block_starts <- c(1, cumsum(block_sizes[-length(block_sizes)]) + 1)
  col_start    <- block_starts[method_index]
  col_end      <- col_start + block_sizes[method_index] - 1
  alpha_samples <- alpha_comm[, col_start:col_end, drop = FALSE]

  if (nsamp > nrow(alpha_samples)) {
    stop("nsamp (", nsamp, ") cannot exceed number of posterior samples (",
         nrow(alpha_samples), ").")
  }
  alpha_samples <- alpha_samples[seq_len(nsamp), , drop = FALSE]
  coef_names    <- colnames(alpha_samples)

  comb_df <- expand.grid(
    variable   = variables_to_scale,
    year_label = class_levels_labels,
    stringsAsFactors = FALSE
  )
  class_label_to_model <- stats::setNames(class_levels_model, class_levels_labels)
  n_comb <- nrow(comb_df)

  predDet <- array(
    NA_real_, dim = c(n_pred, nsamp, n_comb),
    dimnames = list(NULL, NULL,
                    apply(comb_df, 1, function(x) paste(x, collapse = "_")))
  )

  for (c in seq_len(n_comb)) {
    var         <- comb_df$variable[c]
    class_label <- comb_df$year_label[c]
    class_coef  <- class_label_to_model[class_label]
    x           <- scaled_centered_list[[var]]

    for (i in seq_len(nsamp)) {
      eta <- alpha_samples[i, "(Intercept)"]
      if (var %in% coef_names) eta <- eta + alpha_samples[i, var] * x
      if (!is.null(class_coef) && class_coef != "(Intercept)" &&
          class_coef %in% coef_names) {
        eta <- eta + alpha_samples[i, class_coef]
      }
      predDet[, i, c] <- stats::plogis(eta)
    }
  }

  list(predDet = predDet, comb_df = comb_df)
}


#' Run the Full Detection Prediction Pipeline
#'
#' Wrapper that executes the complete detection prediction workflow by
#' sequentially calling \code{\link{build_occ_sequences}},
#' \code{\link{compute_detection}}, and \code{\link{build_occ_summary}}.
#' Mirrors \code{\link{compute_occ_predictions}} for the detection sub-model.
#'
#' @author Thiago Bicudo
#' @seealso \code{\link{compute_occ_predictions}}, \code{\link{compute_detection}},
#'   \code{\link{plot_occ_predictions}}
#' @export
compute_det_predictions <- function(variables_to_scale,
                                    variables_df,
                                    out_bestmodel,
                                    class_levels_labels,
                                    class_levels_model,
                                    nsamp,
                                    n_pred          = 100L,
                                    method_index    = 1L,
                                    variable_labels = NULL) {

  stopifnot(
    "variables_to_scale must be character"  = is.character(variables_to_scale),
    "variables_df must be a data.frame"     = is.data.frame(variables_df),
    "class_levels_labels must be character" = is.character(class_levels_labels),
    "class_levels_model must be character"  = is.character(class_levels_model),
    "class_levels_labels and class_levels_model must have the same length" =
      length(class_levels_labels) == length(class_levels_model),
    "nsamp must be numeric"                 = is.numeric(nsamp),
    "n_pred must be numeric"                = is.numeric(n_pred)
  )

  prep <- build_occ_sequences(
    variables_to_scale = variables_to_scale,
    variables_df       = variables_df,
    n_pred             = n_pred
  )

  pred <- compute_detection(
    out_bestmodel        = out_bestmodel,
    scaled_centered_list = prep$scaled_centered_list,
    variables_to_scale   = variables_to_scale,
    class_levels_labels  = class_levels_labels,
    class_levels_model   = class_levels_model,
    nsamp                = nsamp,
    n_pred               = n_pred,
    method_index         = method_index
  )

  pred_det_df <- build_occ_summary(
    predOcc            = pred$predDet,
    year_levels_labels = class_levels_labels,
    type_levels         = NULL,
    seq_list            = prep$seq_list,
    variable_labels     = variable_labels
  )

  list(
    seq_list             = prep$seq_list,
    scaled_centered_list = prep$scaled_centered_list,
    predDet              = pred$predDet,
    pred_det_df          = pred_det_df
  )
}
