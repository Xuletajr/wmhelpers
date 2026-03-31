#' Run the Full Occupancy Prediction Pipeline
#'
#' Wrapper function that executes the full occupancy prediction workflow by
#' sequentially calling \code{\link{prepare_occ_sequences}},
#' \code{\link{predict_occupancy}}, and
#' \code{\link{summarise_occ_predictions}}.
#'
#' This function standardizes continuous covariates, generates prediction
#' sequences, computes posterior occupancy predictions, and summarizes results
#' into tidy outputs ready for plotting.
#'
#' @param variables_to_scale Character vector of continuous variable names
#'   to be standardized and used for prediction.
#' @param variables_df Data frame containing raw covariate data.
#' @param out_bestmodel Fitted occupancy model object.
#' @param year_levels_model Character or factor vector representing year or
#'   period levels as defined in \code{beta.comm.samples}.
#' @param type_levels Optional character vector for categorical predictors
#'   (e.g., habitat type). Defaults to \code{NULL}.
#' @param nsamp Integer. Number of posterior samples to use.
#' @param n_pred Integer. Number of prediction points per variable.
#'   Defaults to \code{100L}.
#' @param variable_labels Optional named character vector for relabeling
#'   variables in outputs and plots.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{seq_list}}{List of raw sequences for each predictor.}
#'   \item{\code{scaled_centered_list}}{List of scaled/centered sequences.}
#'   \item{\code{predOcc}}{Array of posterior occupancy predictions.}
#'   \item{\code{pred_occ_df}}{Tidy data frame summarizing predictions
#'   (mean and credible intervals).}
#' }
#'
#' @details
#' This function simplifies the occupancy prediction workflow by chaining
#' together three core steps:
#' \enumerate{
#'   \item Sequence preparation and scaling
#'   \item Posterior prediction
#'   \item Summary statistics extraction
#' }
#'
#' It assumes that all required helper functions are available in the package
#' namespace.
#'
#' @author Thiago Bicudo
#'
#' @export
run_occ_predictions <- function(
    variables_to_scale,
    variables_df,
    out_bestmodel,
    year_levels_model,
    type_levels = NULL,
    nsamp,
    n_pred = 100L,
    variable_labels = NULL
) {

  # ── Input validation ───────────────────────────────────────────────
  stopifnot(
    "variables_to_scale must be character" = is.character(variables_to_scale),
    "variables_df must be a data.frame"    = is.data.frame(variables_df),
    "nsamp must be numeric"                = is.numeric(nsamp),
    "n_pred must be numeric"               = is.numeric(n_pred)
  )

  # ── Step 1: Prepare sequences ──────────────────────────────────────
  prep <- prepare_occ_sequences(
    variables_to_scale = variables_to_scale,
    variables_df       = variables_df,
    n_pred             = n_pred
  )

  # ── Step 2: Predict occupancy ──────────────────────────────────────
  pred <- predict_occupancy(
    out_bestmodel         = out_bestmodel,
    scaled_centered_list  = prep$scaled_centered_list,
    year_levels_model     = year_levels_model,
    type_levels           = type_levels,
    nsamp                 = nsamp
  )

  # ── Step 3: Summarise predictions ──────────────────────────────────
  pred_occ_df <- summarise_occ_predictions(
    predOcc         = pred$predOcc,
    seq_list        = prep$seq_list,
    variable_labels = variable_labels
  )

  # ── Output ─────────────────────────────────────────────────────────
  return(list(
    seq_list             = prep$seq_list,
    scaled_centered_list = prep$scaled_centered_list,
    predOcc              = pred$predOcc,
    pred_occ_df          = pred_occ_df
  ))
}
