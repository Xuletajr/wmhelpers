# =============================================================================
# Occupancy prediction pipeline functions
# File: R/msom_occ_predictions.R
# =============================================================================


#' Build Covariate Sequences for Occupancy Prediction
#'
#' Prepares raw and standardised covariate sequences for occupancy prediction
#' by computing the observed range of each predictor variable and returning
#' both the raw sequence (for plot axes) and the scaled sequence (for model
#' prediction).
#'
#' @author Thiago Bicudo
#'
#' @param variables_to_scale Character vector of variable names to process.
#'   Each name must exist as a column in \code{variables_df}.
#' @param variables_df Data frame containing raw covariate data. Must contain
#'   all columns named in \code{variables_to_scale}.
#' @param n_pred Integer. Number of evenly spaced prediction points spanning
#'   the observed range of each variable. Defaults to \code{100L}.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\code{seq_list}}{Named list of raw (unscaled) sequences, one per
#'     variable. Used to populate x-axis values in plots.}
#'   \item{\code{scaled_centered_list}}{Named list of scaled and centred
#'     sequences, one per variable. Passed directly to
#'     \code{\link{compute_occupancy}}.}
#' }
#'
#' @details
#' Scaling uses the mean and standard deviation computed from non-missing
#' values of each variable in \code{variables_df} — the same parameters used
#' to standardise the covariates before model fitting. Variables with zero
#' variance raise an error.
#'
#' @seealso \code{\link{compute_occupancy}}, \code{\link{build_occ_summary}},
#'   \code{\link{compute_occ_predictions}}
#'
#' @export
build_occ_sequences <- function(variables_to_scale,
                                variables_df,
                                n_pred = 100L) {

  # ── Input validation ─────────────────────────────────────────────────────────
  stopifnot(
    "variables_to_scale must be character" = is.character(variables_to_scale),
    "variables_df must be a data.frame"    = is.data.frame(variables_df),
    "n_pred must be numeric"               = is.numeric(n_pred)
  )

  missing_vars <- setdiff(variables_to_scale, names(variables_df))
  if (length(missing_vars) > 0) {
    stop(
      "The following variables were not found in variables_df:\n",
      paste(missing_vars, collapse = ", ")
    )
  }

  seq_list             <- list()
  scaled_centered_list <- list()

  for (var in variables_to_scale) {

    x       <- variables_df[[var]]
    x_no_na <- x[!is.na(x)]

    seq_x <- seq(min(x_no_na), max(x_no_na), length.out = n_pred)

    m    <- mean(x_no_na)
    sd_x <- stats::sd(x_no_na)

    if (sd_x == 0) {
      stop("Variable '", var, "' has zero variance and cannot be scaled.")
    }

    seq_list[[var]]             <- seq_x
    scaled_centered_list[[var]] <- (seq_x - m) / sd_x
  }

  list(
    seq_list             = seq_list,
    scaled_centered_list = scaled_centered_list
  )
}


#' Predict Community-Level Occupancy Across Covariate Gradients
#'
#' Generates posterior occupancy predictions across a sequence of standardised
#' covariate values using community-level beta coefficients from a fitted
#' Multi-Species Occupancy Model (spOccupancy). Supports single- and
#' multi-period models as well as optional categorical type predictors.
#'
#' @author Thiago Bicudo
#'
#' @param out_bestmodel Fitted model object containing
#'   \code{beta.comm.samples} (community-level occupancy coefficient
#'   posterior samples). Typically the output of
#'   \code{spOccupancy::intMsPGOcc()} or \code{spOccupancy::msPGOcc()}.
#' @param scaled_centered_list Named list of scaled covariate sequences,
#'   as returned by \code{\link{build_occ_sequences}}.
#' @param variables_to_scale Character vector of continuous variable names.
#'   Each name must appear in both \code{scaled_centered_list} and the column
#'   names of \code{beta.comm.samples}.
#' @param year_levels_labels Character vector of display labels for year or
#'   survey period levels (e.g. \code{c("2025")} for a single period or
#'   \code{c("2022", "2023")} for a multi-year model). Used in
#'   \code{dimnames} of the output array and in plot titles. Must be parallel
#'   to \code{year_levels_model}.
#' @param year_levels_model Character vector of coefficient names exactly as
#'   they appear in \code{beta.comm.samples} (e.g. \code{c("(Intercept)")}
#'   for a single-period model or \code{c("(Intercept)", "year2023")} for a
#'   multi-year model). The first element is always treated as the reference
#'   level (intercept only, no additional coefficient added). Must be parallel
#'   to \code{year_levels_labels}.
#' @param type_levels Optional character vector of categorical class or habitat
#'   type levels. The first element is treated as the reference level
#'   (intercept only). Defaults to \code{NULL}.
#' @param nsamp Integer. Number of posterior samples to use for prediction.
#'   Must not exceed \code{nrow(out_bestmodel$beta.comm.samples)}.
#' @param n_pred Integer. Number of prediction points per variable.
#'   Defaults to \code{100L}.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\code{predOcc}}{3-dimensional array \code{[n_pred x nsamp x
#'     n_combinations]} of posterior occupancy probabilities on the
#'     probability scale (after \code{plogis}). The third dimension has named
#'     \code{dimnames} built from \code{variable_year_label} (and
#'     \code{_type} if applicable).}
#'   \item{\code{comb_df}}{Data frame of all variable x year_label (x type)
#'     combinations used to build the array.}
#' }
#'
#' @details
#' \code{year_levels_labels} and \code{year_levels_model} must be parallel
#' vectors of equal length. Labels are used for \code{dimnames} and plotting;
#' model names are used only to look up coefficients in
#' \code{beta.comm.samples}.
#'
#' For single-period models, use \code{year_levels_labels = c("2025")} and
#' \code{year_levels_model = c("(Intercept)")}. The year label will be
#' suppressed from output titles by \code{\link{build_occ_summary}} when only
#' one period exists.
#'
#' @seealso \code{\link{build_occ_sequences}}, \code{\link{build_occ_summary}},
#'   \code{\link{compute_occ_predictions}}
#'
#' @importFrom stats plogis
#'
#' @export
compute_occupancy <- function(out_bestmodel,
                              scaled_centered_list,
                              variables_to_scale,
                              year_levels_labels,
                              year_levels_model,
                              type_levels = NULL,
                              nsamp,
                              n_pred      = 100L) {

  # ── Input validation ─────────────────────────────────────────────────────────
  stopifnot(
    "out_bestmodel must contain beta.comm.samples" =
      !is.null(out_bestmodel$beta.comm.samples),
    "scaled_centered_list must be a list"          = is.list(scaled_centered_list),
    "variables_to_scale must be character"         = is.character(variables_to_scale),
    "year_levels_labels must be character"         = is.character(year_levels_labels),
    "year_levels_model must be character"          = is.character(year_levels_model),
    "year_levels_labels and year_levels_model must have the same length" =
      length(year_levels_labels) == length(year_levels_model),
    "nsamp must be numeric"                        = is.numeric(nsamp),
    "n_pred must be numeric"                       = is.numeric(n_pred)
  )

  beta_samples <- out_bestmodel$beta.comm.samples

  if (nsamp > nrow(beta_samples)) {
    stop(
      "nsamp (", nsamp, ") cannot exceed number of posterior samples (",
      nrow(beta_samples), ")."
    )
  }

  beta_samples <- beta_samples[seq_len(nsamp), , drop = FALSE]
  coef_names   <- colnames(beta_samples)

  # ── Build combinations using display labels ───────────────────────────────────
  if (is.null(type_levels)) {
    comb_df <- expand.grid(
      variable   = variables_to_scale,
      year_label = year_levels_labels,
      stringsAsFactors = FALSE
    )
  } else {
    comb_df <- expand.grid(
      variable   = variables_to_scale,
      year_label = year_levels_labels,
      type       = type_levels,
      stringsAsFactors = FALSE
    )
  }

  year_label_to_model <- stats::setNames(year_levels_model, year_levels_labels)
  n_comb              <- nrow(comb_df)

  predOcc <- array(
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

    var        <- comb_df$variable[c]
    year_label <- comb_df$year_label[c]
    year_coef  <- year_label_to_model[year_label]
    x          <- scaled_centered_list[[var]]
    type       <- if (!is.null(type_levels)) comb_df$type[c] else NULL

    for (i in seq_len(nsamp)) {

      eta <- beta_samples[i, "(Intercept)"]

      # Continuous variable effect
      if (var %in% coef_names) {
        eta <- eta + beta_samples[i, var] * x
      }

      # Year effect — skip if reference level (intercept only)
      if (!is.null(year_coef) &&
          year_coef != "(Intercept)" &&
          year_coef %in% coef_names) {
        eta <- eta + beta_samples[i, year_coef]
      }

      # Type effect — skip if reference level
      if (!is.null(type)) {
        if (type != type_levels[1] && type %in% coef_names) {
          eta <- eta + beta_samples[i, type]
        }
      }

      predOcc[, i, c] <- stats::plogis(eta)
    }
  }

  list(predOcc = predOcc, comb_df = comb_df)
}


#' Summarise Posterior Occupancy Predictions into a Tidy Data Frame
#'
#' Converts the 3-dimensional prediction array returned by
#' \code{\link{compute_occupancy}} into a tidy \code{data.frame} with
#' posterior mean, standard deviation, and 95\% credible intervals for each
#' predictor x year (x class) combination.
#'
#' @author Thiago Bicudo
#'
#' @param predOcc 3-dimensional array \code{[n_pred x nsamp x n_combinations]}
#'   as returned by \code{\link{compute_occupancy}}. The third dimension must
#'   have named \code{dimnames} built from display labels.
#' @param year_levels_labels Character vector of display labels used for
#'   year/period levels. Used to parse combination names and suppressed from
#'   output when only one level exists (single-period models).
#' @param type_levels Optional character vector of class/type display labels.
#'   Defaults to \code{NULL}.
#' @param seq_list Named list of raw (unscaled) covariate sequences as
#'   returned by \code{\link{build_occ_sequences}}. Used to populate the
#'   \code{seq_data} column (x-axis values for plots). Defaults to
#'   \code{NULL}, in which case integer indices are used instead.
#' @param variable_labels Optional named character vector mapping short
#'   variable names to human-readable labels for plot axes
#'   (e.g. \code{c(twi = "Topographic Wetness Index")}). Defaults to
#'   \code{NULL}, in which case the short variable name is used as-is.
#'
#' @return A \code{data.frame} with one row per prediction point per
#'   combination, and columns:
#' \describe{
#'   \item{\code{variable}}{Short variable name.}
#'   \item{\code{expanded_variable}}{Human-readable label from
#'     \code{variable_labels}, or \code{variable} if not provided.}
#'   \item{\code{year}}{ Year/period display label, or \code{NA} for
#'     single-period models.}
#'   \item{\code{class}}{Type/class label, or \code{NA} if
#'     \code{type_levels} is \code{NULL}.}
#'   \item{\code{mean}}{Posterior mean occupancy probability.}
#'   \item{\code{sd}}{Posterior standard deviation.}
#'   \item{\code{ciL}}{Lower bound of 95\% credible interval (2.5th
#'     percentile).}
#'   \item{\code{ciU}}{Upper bound of 95\% credible interval (97.5th
#'     percentile).}
#'   \item{\code{seq_data}}{Raw covariate values for the x-axis.}
#' }
#'
#' @details
#' Column names are aligned with the expected input of
#' \code{\link{plot_occ_predictions}}: \code{seq_data} for the x-axis and
#' \code{ciL}/\code{ciU} for credible interval bounds.
#'
#' Variable names in \code{predOcc} dimnames may contain underscores (e.g.
#' \code{"dist_built_2025"}). The function matches year and type labels from
#' the end of the combination name string to correctly isolate the variable
#' name prefix.
#'
#' @seealso \code{\link{build_occ_sequences}}, \code{\link{compute_occupancy}},
#'   \code{\link{compute_occ_predictions}}, \code{\link{plot_occ_predictions}}
#'
#' @importFrom stats quantile sd
#'
#' @export
build_occ_summary <- function(predOcc,
                              year_levels_labels,
                              type_levels     = NULL,
                              seq_list        = NULL,
                              variable_labels = NULL) {

  # ── Input validation ─────────────────────────────────────────────────────────
  stopifnot(
    "predOcc must be an array"             = is.array(predOcc),
    "predOcc must be 3-dimensional"        = length(dim(predOcc)) == 3,
    "year_levels_labels must be character" = is.character(year_levels_labels)
  )

  n_pred     <- dim(predOcc)[1]
  comb_names <- dimnames(predOcc)[[3]]

  if (is.null(comb_names)) {
    stop("predOcc must have named third dimension (combination labels).")
  }

  has_type    <- !is.null(type_levels)
  single_year <- length(year_levels_labels) == 1

  # ── Helper: parse combination name into variable / year / type ───────────────
  # Variable names may contain underscores — match year and type from end.
  parse_comb_name <- function(comb_name) {

    if (has_type) {
      matched_type <- type_levels[
        sapply(type_levels, function(t) endsWith(comb_name, paste0("_", t)))
      ]
      type <- if (length(matched_type) > 0) matched_type[1] else NA_character_

      matched_year <- year_levels_labels[
        sapply(year_levels_labels, function(y)
          grepl(paste0("_", y, "_"), comb_name, fixed = TRUE) |
            grepl(paste0("_", y, "$"), comb_name))
      ]
      year <- if (length(matched_year) > 0) matched_year[1] else NA_character_

      variable <- sub(paste0("_", year, "_", type, "$"), "", comb_name)

    } else {
      matched_year <- year_levels_labels[
        sapply(year_levels_labels, function(y)
          endsWith(comb_name, paste0("_", y)))
      ]
      year     <- if (length(matched_year) > 0) matched_year[1] else NA_character_
      variable <- sub(paste0("_", year, "$"), "", comb_name)
      type     <- NA_character_
    }

    display_year <- if (single_year) NA_character_ else year

    list(variable = variable, year = display_year, type = type)
  }

  # ── Build output ──────────────────────────────────────────────────────────────
  out_list <- vector("list", length(comb_names))

  for (i in seq_along(comb_names)) {

    parsed   <- parse_comb_name(comb_names[i])
    variable <- parsed$variable
    year     <- parsed$year
    type     <- parsed$type

    mat <- predOcc[, , i]

    mean_pred <- apply(mat, 1, mean)
    sd_pred   <- apply(mat, 1, stats::sd)
    ciL       <- apply(mat, 1, stats::quantile, probs = 0.025)
    ciU       <- apply(mat, 1, stats::quantile, probs = 0.975)

    seq_data <- if (!is.null(seq_list) && variable %in% names(seq_list)) {
      seq_list[[variable]]
    } else {
      seq_len(n_pred)
    }

    expanded_variable <- if (!is.null(variable_labels) &&
                             variable %in% names(variable_labels)) {
      variable_labels[[variable]]
    } else {
      variable
    }

    out_list[[i]] <- data.frame(
      variable          = rep(variable, n_pred),
      expanded_variable = rep(expanded_variable, n_pred),
      year              = rep(year, n_pred),
      class             = rep(type, n_pred),
      mean              = mean_pred,
      sd                = sd_pred,
      ciL               = ciL,
      ciU               = ciU,
      seq_data          = seq_data,
      stringsAsFactors  = FALSE
    )
  }

  out_df           <- do.call(rbind, out_list)
  rownames(out_df) <- NULL

  out_df
}


#' Run the Full Occupancy Prediction Pipeline
#'
#' Wrapper that executes the complete occupancy prediction workflow by
#' sequentially calling \code{\link{build_occ_sequences}},
#' \code{\link{compute_occupancy}}, and \code{\link{build_occ_summary}}.
#'
#' @author Thiago Bicudo
#'
#' @param variables_to_scale Character vector of continuous variable names
#'   to be standardised and used for prediction. Each name must exist as a
#'   column in \code{variables_df} and in \code{beta.comm.samples}.
#' @param variables_df Data frame containing raw covariate data.
#' @param out_bestmodel Fitted occupancy model object containing
#'   \code{beta.comm.samples}.
#' @param year_levels_labels Character vector of display labels for year or
#'   survey period levels (e.g. \code{c("2025")}). Must be parallel to
#'   \code{year_levels_model}.
#' @param year_levels_model Character vector of coefficient names exactly as
#'   they appear in \code{beta.comm.samples} (e.g. \code{c("(Intercept)")}
#'   for single-period models). Must be parallel to \code{year_levels_labels}.
#' @param type_levels Optional character vector for categorical predictors
#'   (e.g., habitat types). The first element is treated as the reference
#'   level. Defaults to \code{NULL}.
#' @param nsamp Integer. Number of posterior samples to use for prediction.
#' @param n_pred Integer. Number of prediction points per variable.
#'   Defaults to \code{100L}.
#' @param variable_labels Optional named character vector mapping short
#'   variable names to human-readable axis labels
#'   (e.g. \code{c(twi = "Topographic Wetness Index")}).
#'   Defaults to \code{NULL}.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{seq_list}}{Named list of raw covariate sequences per
#'     variable.}
#'   \item{\code{scaled_centered_list}}{Named list of scaled and centred
#'     sequences per variable.}
#'   \item{\code{predOcc}}{3-dimensional array of posterior occupancy
#'     probabilities \code{[n_pred x nsamp x n_combinations]}.}
#'   \item{\code{pred_occ_df}}{Tidy data frame of posterior summaries
#'     (mean, SD, credible intervals) ready for \code{\link{plot_occ_predictions}}.}
#' }
#'
#' @details
#' \code{year_levels_labels} and \code{year_levels_model} must be parallel
#' vectors of equal length. Labels control \code{dimnames} and plot display;
#' model names are used only to look up coefficients in
#' \code{beta.comm.samples}.
#'
#' For single-period models use \code{year_levels_labels = c("2025")} and
#' \code{year_levels_model = c("(Intercept)")}. Year labels will be
#' suppressed automatically from plot titles.
#'
#' Typical usage:
#' \preformatted{
#' results <- compute_occ_predictions(
#'   variables_to_scale = c("twi", "dist_built"),
#'   variables_df       = variables,
#'   out_bestmodel      = out_bestmodel,
#'   year_levels_labels = c("2025"),
#'   year_levels_model  = c("(Intercept)"),
#'   type_levels        = NULL,
#'   nsamp              = nsamp,
#'   n_pred             = 100,
#'   variable_labels    = c(
#'     twi        = "Topographic Wetness Index",
#'     dist_built = "Distance to development (buildings)"
#'   )
#' )
#'
#' plot_list <- plot_occ_predictions(pred_occ_df = results$pred_occ_df)
#' }
#'
#' @seealso \code{\link{build_occ_sequences}}, \code{\link{compute_occupancy}},
#'   \code{\link{build_occ_summary}}, \code{\link{plot_occ_predictions}}
#'
#' @export
compute_occ_predictions <- function(variables_to_scale,
                                    variables_df,
                                    out_bestmodel,
                                    year_levels_labels,
                                    year_levels_model,
                                    type_levels     = NULL,
                                    nsamp,
                                    n_pred          = 100L,
                                    variable_labels = NULL) {

  # ── Input validation ─────────────────────────────────────────────────────────
  stopifnot(
    "variables_to_scale must be character" = is.character(variables_to_scale),
    "variables_df must be a data.frame"    = is.data.frame(variables_df),
    "year_levels_labels must be character" = is.character(year_levels_labels),
    "year_levels_model must be character"  = is.character(year_levels_model),
    "year_levels_labels and year_levels_model must have the same length" =
      length(year_levels_labels) == length(year_levels_model),
    "nsamp must be numeric"                = is.numeric(nsamp),
    "n_pred must be numeric"               = is.numeric(n_pred)
  )

  # ── Step 1: Build covariate sequences ────────────────────────────────────────
  prep <- build_occ_sequences(
    variables_to_scale = variables_to_scale,
    variables_df       = variables_df,
    n_pred             = n_pred
  )

  # ── Step 2: Compute posterior occupancy predictions ───────────────────────────
  pred <- compute_occupancy(
    out_bestmodel        = out_bestmodel,
    scaled_centered_list = prep$scaled_centered_list,
    variables_to_scale   = variables_to_scale,
    year_levels_labels   = year_levels_labels,
    year_levels_model    = year_levels_model,
    type_levels          = type_levels,
    nsamp                = nsamp,
    n_pred               = n_pred
  )

  # ── Step 3: Summarise into tidy data frame ────────────────────────────────────
  pred_occ_df <- build_occ_summary(
    predOcc            = pred$predOcc,
    year_levels_labels = year_levels_labels,
    type_levels        = type_levels,
    seq_list           = prep$seq_list,
    variable_labels    = variable_labels
  )

  list(
    seq_list             = prep$seq_list,
    scaled_centered_list = prep$scaled_centered_list,
    predOcc              = pred$predOcc,
    pred_occ_df          = pred_occ_df
  )
}
