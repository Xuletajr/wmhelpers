# =============================================================================
# wmhelper — MSOM Inference & Plotting Functions
#
# Functions in this file (in pipeline order):
#
#  DIAGNOSTICS
#   1. check_rhat()              - R-hat per parameter, returns problematic table
#   2. plot_rhat_hist()          - histogram of R-hat values (beta + alpha)
#   3. check_ess()               - ESS per parameter, returns problematic table
#   4. plot_ess_hist()           - histogram of ESS values (beta + alpha)
#
#  PARAMETER SUMMARIES
#   5. extract_prob()            - occupancy OR detection probability per species
#   6. build_param_table()       - beta or alpha species-level summary table
#   7. get_significant_vars()    - variables whose 95% CRI excludes zero
#
#  COMPARISON PLOT
#   8. plot_occ_vs_det()         - scatter of mean occupancy vs mean detection
#
#  DETECTION PREDICTION
#   9. predict_detection_julian() - community-level detection over Julian days
#  10. plot_detection_julian()    - plot the detection curve
#
#  OCCUPANCY PREDICTION PIPELINE (continuous covariates)
#  11. prepare_occ_sequences()   - raw + scaled sequences per variable
#  12. predict_occupancy()       - fill 3-D prediction array from MCMC
#  13. summarise_occ_predictions()- tidy stats data frame
#  14. run_occ_predictions()     - master wrapper for 11-13
#
#  OCCUPANCY PLOTS
#  15. plot_occ_predictions()    - list of ggplot objects
#  16. export_occ_plots()        - individual PNGs + multiplot
#
#  PSI (SITE-LEVEL OCCUPANCY)
#  17. extract_psi_stats()       - mean/sd/CRI per species x site
#  18. export_psi_tables()       - write CSVs for each statistic
#
# =============================================================================
# WHAT TO KEEP AS PLAIN SCRIPT (not wrapped in functions):
#   - setwd / rm / gc / library calls
#   - yaml::yaml.load_file()  and  load() / readRDS()
#   - out_bestmodel$x.names / dim() / identical() inspection calls
#   - sink() model summary export  (project-specific boilerplate)
#   - variables rename block  (project-specific column mapping)
#   - year_levels / type_levels / variable_labels definitions
#   - print(plot_list[[i]])  interactive checks
# =============================================================================


# ── Dependencies (declared here for package DESCRIPTION Imports:) ─────────────
# dplyr, tidyr, tibble, purrr, stringr, ggplot2, ggpubr, rlang


# =============================================================================
# 1. check_rhat
# =============================================================================

#' Check R-hat Convergence Diagnostic
#'
#' Evaluates R-hat (Gelman-Rubin) values for a set of MCMC parameters and
#' returns a table of any parameters exceeding the acceptable threshold.
#' Issues a \code{message()} when all values are acceptable and a
#' \code{warning()} when any exceed the threshold.
#'
#' @param rhat_values Named numeric vector of R-hat values
#'   (e.g. \code{out_bestmodel$rhat$beta}).
#' @param samples Matrix of MCMC samples whose column names identify
#'   parameters (e.g. \code{out_bestmodel$beta.samples}). Used only to
#'   recover parameter names.
#' @param parameter_name Character label used in messages
#'   (e.g. \code{"BETA"} or \code{"ALPHA"}).
#' @param threshold Numeric. R-hat values at or above this value are flagged.
#'   Defaults to \code{1.1}.
#'
#' @return A \code{data.frame} with columns \code{Parameter} and \code{Rhat}
#'   for all problematic parameters, or \code{NULL} if none are found.
#'
#' @export
check_rhat <- function(rhat_values,
                       samples,
                       parameter_name,
                       threshold = 1.1) {

  stopifnot(
    "rhat_values must be a numeric vector"  = is.numeric(rhat_values),
    "samples must have column names"        = !is.null(colnames(samples)),
    "parameter_name must be a string"       = is.character(parameter_name)
  )

  problematic_idx <- which(rhat_values >= threshold)

  if (length(problematic_idx) > 0) {

    warning(
      "MCMC chains have not mixed well for the following ",
      parameter_name, " parameters. Consider refitting the model.",
      call. = FALSE
    )

    result <- data.frame(
      Parameter = colnames(samples)[problematic_idx],
      Rhat      = round(rhat_values[problematic_idx], 3)
    )
    print(result)
    return(result)

  } else {
    message(
      "R-hat acceptable for ", parameter_name, " parameters ",
      "(max = ", round(max(rhat_values, na.rm = TRUE), 4), ")."
    )
    return(NULL)
  }
}


# =============================================================================
# 2. plot_rhat_hist
# =============================================================================

#' Plot Histogram of R-hat Values
#'
#' Produces a faceted histogram of R-hat values for beta (occupancy) and
#' alpha (detection) parameters and optionally exports it to PNG.
#'
#' @param out_bestmodel Fitted model object with \code{rhat$beta} and
#'   \code{rhat$alpha} numeric vectors.
#' @param output_path Character. Full path including filename for PNG export.
#'   If \code{NULL} (default), the plot is printed to the active device only.
#' @param color_fill Character. Fill and border colour for histogram bars.
#'   Defaults to \code{"#CCCCCC"}.
#' @param color_text Character. Axis text colour. Defaults to \code{"#1B3A5C"}.
#' @param color_grid Character. Panel grid colour. Defaults to \code{"#E8EEF4"}.
#' @param width,height,dpi Numeric. PNG dimensions in inches and resolution.
#'
#' @return Invisibly returns the \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot geom_histogram facet_wrap labs theme_bw theme
#'   element_text element_rect element_line aes
#' @importFrom tibble tibble
#'
#' @export
plot_rhat_hist <- function(out_bestmodel,
                           output_path  = NULL,
                           color_fill   = "#CCCCCC",
                           color_text   = "#1B3A5C",
                           color_grid   = "#E8EEF4",
                           width        = 12,
                           height       = 9,
                           dpi          = 300) {

  p <- tibble::tibble(
    rhat  = c(out_bestmodel$rhat$beta,  out_bestmodel$rhat$alpha),
    param = c(rep("beta",  length(out_bestmodel$rhat$beta)),
              rep("alpha", length(out_bestmodel$rhat$alpha)))
  ) %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(
      ggplot2::aes(x = rhat),
      color = color_fill, fill = color_fill
    ) +
    ggplot2::facet_wrap(~ param, ncol = 2, scales = "fixed") +
    ggplot2::labs(x = "Rhat (Gelman-Rubin diagnostic)") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title   = ggplot2::element_text(size = 18, color = color_text),
      axis.text    = ggplot2::element_text(size = 14, color = color_text),
      strip.text.x = ggplot2::element_text(size = 16, color = color_text,
                                           face = "bold"),
      strip.background = ggplot2::element_rect(fill = color_grid),
      panel.grid   = ggplot2::element_line(color = color_grid)
    )

  if (!is.null(output_path)) {
    png(output_path, width = width, height = height,
        units = "in", res = dpi, bg = "white")
    print(p)
    dev.off()
    message("R-hat histogram saved: ", basename(output_path))
  } else {
    print(p)
  }
return(p)
}


# =============================================================================
# 3. check_ess
# =============================================================================

#' Check Effective Sample Size (ESS)
#'
#' Identifies parameters with ESS below a minimum threshold and optionally
#' exports a CSV of problematic parameters.
#'
#' @param ess_values Named numeric vector of ESS values
#'   (e.g. \code{out_bestmodel$ESS$beta}).
#' @param parameter_name Character label for messages (e.g. \code{"BETA"}).
#' @param threshold Numeric. ESS values below this are flagged.
#'   Defaults to \code{450} (150 per chain × 3 chains).
#' @param output_path Character. If provided, the problematic table is written
#'   to this CSV path. Defaults to \code{NULL}.
#'
#' @return A \code{data.frame} with columns \code{Parameter} and \code{ESS}
#'   for all flagged parameters, or \code{NULL} if none found.
#'
#' @export
check_ess <- function(ess_values,
                      parameter_name,
                      threshold   = 450,
                      output_path = NULL) {

  stopifnot(
    "ess_values must be a named numeric vector" =
      is.numeric(ess_values) && !is.null(names(ess_values))
  )

  problematic <- sort(ess_values[ess_values < threshold])

  if (length(problematic) == 0) {
    message(
      "All ESS values for ", parameter_name,
      " are above the minimum threshold (", threshold, ")."
    )
    return(NULL)
  }

  result <- data.frame(
    Parameter = names(problematic),
    ESS       = round(as.numeric(problematic))
  )

  message(
    length(problematic), " ", parameter_name,
    " parameter(s) with ESS < ", threshold, "."
  )

  if (!is.null(output_path)) {
    readr::write_csv(result, output_path)
    message("Problematic ESS table saved: ", basename(output_path))
  }

  result
}


# =============================================================================
# 4. plot_ess_hist
# =============================================================================

#' Plot Histogram of Effective Sample Size (ESS) Values
#'
#' Produces a faceted histogram of ESS values for beta and alpha parameters
#' and optionally exports to PNG.
#'
#' @inheritParams plot_rhat_hist
#' @param fig_font Character. Font family for axis text. Defaults to \code{""}
#'
#' @return Invisibly returns the \code{ggplot} object.
#'
#' @export
plot_ess_hist <- function(out_bestmodel,
                          output_path = NULL,
                          color_fill  = "#CCCCCC",
                          color_text  = "#1B3A5C",
                          color_grid  = "#E8EEF4",
                          fig_font    = "",
                          width       = 12,
                          height      = 9,
                          dpi         = 300) {

  p <- tibble::tibble(
    ess   = c(out_bestmodel$ESS$beta,  out_bestmodel$ESS$alpha),
    param = c(rep("beta",  length(out_bestmodel$ESS$beta)),
              rep("alpha", length(out_bestmodel$ESS$alpha)))
  ) %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(
      ggplot2::aes(x = ess),
      color = color_fill, fill = color_fill
    ) +
    ggplot2::facet_wrap(~ param, ncol = 2, scales = "fixed") +
    ggplot2::labs(x = "Effective Sample Size (ESS)") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title   = ggplot2::element_text(size = 22, color = color_text,
                                           family = fig_font),
      axis.text.x  = ggplot2::element_text(size = 18, color = color_text,
                                           family = fig_font,
                                           angle = 45, hjust = 1),
      axis.text.y  = ggplot2::element_text(size = 18, color = color_text,
                                           family = fig_font),
      strip.text.x = ggplot2::element_text(size = 16, color = color_text,
                                           family = fig_font, face = "bold"),
      strip.background = ggplot2::element_rect(fill = color_grid),
      panel.grid   = ggplot2::element_line(color = color_grid)
    )

  if (!is.null(output_path)) {
    png(output_path, width = width, height = height,
        units = "in", res = dpi, bg = "white")
    print(p)
    dev.off()
    message("ESS histogram saved: ", basename(output_path))
  } else {
    print(p)
  }
return(p)
}


# =============================================================================
# 5. extract_prob
# =============================================================================

#' Extract Species-Level Occupancy or Detection Probabilities
#'
#' Back-transforms MCMC samples for one or more categorical classes to the
#' probability scale via \code{plogis}, then summarises across samples
#' (mean, SD, 95\% CRI). Works for both occupancy (beta) and detection
#' (alpha) parameters.
#'
#' @param samples Matrix of MCMC samples (\code{[nsamp × parameters]}).
#'   For occupancy use \code{out_bestmodel$beta.samples}; for detection use
#'   \code{out_bestmodel$alpha.samples}.
#' @param class_list Named list mapping class labels to their additional
#'   beta/alpha coefficient name, or \code{NULL} for intercept-only.
#'   Example (intercept only):
#'   \code{list(survey2025 = NULL)}
#'   Example (two years, 2022 as reference):
#'   \code{list(survey2022 = NULL, survey2023 = "year2023")}
#' @param spp_reference A data frame with at least a \code{species} column
#'   used to left-join and validate species names
#'   (e.g. \code{spp_list_ea}).
#'
#' @return A tidy \code{tibble} with columns \code{species},
#'   \code{param_spp}, \code{year}, \code{class}, \code{mean},
#'   \code{sd}, \code{ciL}, \code{ciU}.
#'
#' @details
#' The function always uses the \code{(Intercept)} column from \code{samples}
#' as the base. For each class in \code{class_list}, if the value is not
#' \code{NULL}, the corresponding coefficient column is added before applying
#' \code{plogis}.
#'
#' @importFrom dplyr bind_rows group_by summarise mutate left_join select
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_remove str_extract
#' @importFrom purrr imap
#'
#' @export
extract_prob <- function(samples,
                         class_list,
                         spp_reference) {

  stopifnot(
    "samples must be a matrix with column names" =
      is.matrix(samples) && !is.null(colnames(samples)),
    "class_list must be a named list" =
      is.list(class_list) && !is.null(names(class_list)),
    "spp_reference must be a data frame with a 'species' column" =
      is.data.frame(spp_reference) && "species" %in% names(spp_reference)
  )

  long_list <- purrr::imap(class_list, function(extra_coef, class_label) {

    # Base: intercept
    base <- samples[, grep("Intercept", colnames(samples)), drop = FALSE]

    # Add class-specific coefficient if provided
    if (!is.null(extra_coef)) {
      if (!extra_coef %in% colnames(samples)) {
        stop("Column '", extra_coef, "' not found in samples.", call. = FALSE)
      }
      base <- base + samples[, extra_coef, drop = FALSE]
    }

    # Back-transform and reshape
    plogis(base) %>%
      tibble::as_tibble() %>%
      tidyr::pivot_longer(
        cols      = dplyr::everything(),
        names_to  = "param_spp",
        values_to = "estimates"
      ) %>%
      dplyr::mutate(class = class_label)
  })

  dplyr::bind_rows(long_list) %>%
    dplyr::group_by(class, param_spp) %>%
    dplyr::summarise(
      mean = mean(estimates),
      sd   = sd(estimates),
      ciL  = quantile(estimates, probs = 0.025),
      ciU  = quantile(estimates, probs = 0.975),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      species = stringr::str_remove(param_spp, "^\\(Intercept\\)-"),
      year    = stringr::str_extract(class, "\\d{4}"),
      .after  = param_spp
    ) %>%
    dplyr::left_join(
      spp_reference %>% dplyr::select(species),
      by = "species"
    ) %>%
    dplyr::select(species, param_spp, year, class, dplyr::everything())
}


# =============================================================================
# 6. build_param_table
# =============================================================================

#' Build Species-Level Parameter Summary Table
#'
#' Creates a tidy summary table (mean, SD, 95\% CRI, effect direction) for
#' either beta (occupancy) or alpha (detection) species-level MCMC samples.
#' Identical structure for both — select the right \code{samples} matrix.
#'
#' @param samples Matrix of MCMC samples (e.g.
#'   \code{out_bestmodel$beta.samples} or
#'   \code{out_bestmodel$alpha.samples}).
#' @param sp_names Character vector of species names in the same order as
#'   \code{out_bestmodel$sp.names}.
#' @param n_covariates Integer. Number of covariates in the design matrix
#'   (\code{ncol(out_bestmodel$X)} for beta;
#'   \code{ncol(out_bestmodel$X.p)} for alpha).
#' @param output_path Character. If provided, the table is exported to this
#'   CSV path. Defaults to \code{NULL}.
#'
#' @return A \code{data.frame} with columns \code{param_spp_code},
#'   \code{species}, \code{Mean}, \code{SD}, \code{CRI_low},
#'   \code{CRI_high}, \code{Effect}.
#'
#' @details
#' The \code{Effect} column is set to \code{""} for intercept rows
#' (first \code{length(sp_names)} rows) and to \code{"Positive"},
#' \code{"Negative"}, or \code{"No effect"} for covariate rows.
#'
#' @export
build_param_table <- function(samples,
                              sp_names,
                              n_covariates,
                              output_path = NULL) {

  stopifnot(
    "samples must be a matrix"           = is.matrix(samples),
    "sp_names must be a character vector" = is.character(sp_names),
    "n_covariates must be a positive integer" =
      is.numeric(n_covariates) && n_covariates > 0
  )

  smry   <- summary(samples)
  stats  <- round(smry$statistics, 3)
  quants <- round(smry$quantiles,  3)

  result <- cbind(
    param_spp_code = rownames(stats),
    species        = rep(sp_names, n_covariates),
    stats[,  c("Mean", "SD")],
    quants[, c("2.5%", "97.5%")]
  ) %>%
    as.data.frame() %>%
    dplyr::rename(CRI_low = `2.5%`, CRI_high = `97.5%`) %>%
    dplyr::mutate(
      CRI_low  = as.numeric(CRI_low),
      CRI_high = as.numeric(CRI_high),
      Effect   = dplyr::case_when(
        CRI_low > 0 & CRI_high > 0 ~ "Positive",
        CRI_low < 0 & CRI_high < 0 ~ "Negative",
        TRUE                        ~ "No effect"
      )
    )

  # Blank out Effect for intercept rows
  result[seq_len(length(sp_names)), "Effect"] <- ""

  if (!is.null(output_path)) {
    readr::write_csv(result, output_path)
    message("Parameter table saved: ", basename(output_path))
  }

  result
}


# =============================================================================
# 7. get_significant_vars
# =============================================================================

#' Identify Significant Variables from Community-Level Beta Samples
#'
#' Returns variable names whose 95\% credible interval (across community-level
#' beta samples) entirely excludes zero. Used to select variables for
#' multiplots and downstream highlighting.
#'
#' @param beta_comm_samples Matrix of community-level beta MCMC samples
#'   (\code{out_bestmodel$beta.comm.samples}).
#' @param probs Numeric vector of length 2 specifying the lower and upper
#'   quantile bounds. Defaults to \code{c(0.025, 0.975)}.
#' @param exclude Character vector of column names to always exclude
#'   (e.g. intercept). Defaults to \code{"(Intercept)"}.
#'
#' @return Character vector of variable names with CRI excluding zero.
#'
#' @export
get_significant_vars <- function(beta_comm_samples,
                                 probs   = c(0.025, 0.975),
                                 exclude = "(Intercept)") {

  stopifnot(
    "beta_comm_samples must be a matrix" = is.matrix(beta_comm_samples)
  )

  t(apply(beta_comm_samples, 2, quantile, probs = probs)) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("variable") %>%
    dplyr::rename(lower = 2, upper = 3) %>%
    dplyr::filter(!variable %in% exclude) %>%
    dplyr::filter(
      (lower > 0 & upper > 0) |
      (lower < 0 & upper < 0)
    ) %>%
    dplyr::pull(variable)
}


# =============================================================================
# 8. plot_occ_vs_det
# =============================================================================

#' Plot Mean Occupancy vs Mean Detection Probability
#'
#' Scatter plot of species-level mean occupancy (x) against mean detection (y)
#' with a linear regression line, R², and p-value annotation.
#'
#' @param export_occ_prob Tidy data frame from \code{\link{extract_prob}} for
#'   occupancy.
#' @param export_det_prob Tidy data frame from \code{\link{extract_prob}} for
#'   detection.
#' @param color_line Character. Colour of the regression line.
#'   Defaults to \code{"#1B3A5C"}.
#' @param output_path Character. If provided, exported as PNG.
#'   Defaults to \code{NULL}.
#' @param width,height,dpi Numeric. PNG dimensions in inches and resolution.
#'
#' @return Invisibly returns the fitted \code{lm} object.
#'
#' @export
plot_occ_vs_det <- function(export_occ_prob,
                            export_det_prob,
                            color_line  = "#1B3A5C",
                            output_path = NULL,
                            width       = 9,
                            height      = 9,
                            dpi         = 600) {

  x_vals <- export_occ_prob %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(mean = mean(mean), .groups = "drop") %>%
    dplyr::pull(mean)

  y_vals <- export_det_prob %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(mean = mean(mean), .groups = "drop") %>%
    dplyr::pull(mean)

  fit <- lm(y_vals ~ x_vals)
  r2  <- round(summary(fit)$r.squared, 3)
  pv  <- format.pval(summary(fit)$coefficients[2, 4], digits = 3, eps = 0.001)

  # Internal draw function — called once per device so both screen
  # and file get identical output without code duplication
  draw <- function() {
    plot(
      x_vals, y_vals,
      pch  = 16, xlim = c(0, 1), ylim = c(0, 1),
      xlab = "Occurrence probability",
      ylab = "Detection probability"
    )
    abline(fit, lwd = 2, col = color_line)
    text(0.1, 0.90, paste0("R² = ", r2), pos = 4, cex = 1.2)
    text(0.1, 0.85, paste0("p = ",  pv), pos = 4, cex = 1.2)
  }

  # Always draw to screen first
  draw()

  # Then also save to PNG if a path was provided
  if (!is.null(output_path)) {
    png(output_path, width = width, height = height,
        units = "in", res = dpi, pointsize = 12)
    draw()
    dev.off()
    message("Occ vs det plot saved: ", basename(output_path))
  }

  invisible(fit)
}




# =============================================================================
# 9. predict_detection_julian
# =============================================================================

#' Predict Community-Level Detection Probability Across Julian Days
#'
#' Fills a \code{[n_pred × nsamp]} array with back-transformed detection
#' probabilities using a linear + quadratic Julian day effect from community-
#' level alpha samples.
#'
#' @param out_bestmodel Fitted model object.
#' @param julian_day Matrix of Julian day values (sites × surveys). Raw
#'   values used to define the prediction range.
#' @param m_julian Numeric. Mean Julian day used for centering.
#' @param sd_julian Numeric. SD of Julian days used for scaling.
#' @param linear_coef Character. Name of the linear Julian day coefficient
#'   in \code{alpha.comm.samples}. Defaults to \code{"julian_std"}.
#' @param quadratic_coef Character or \code{NULL}. Name of the quadratic
#'   coefficient. Pass \code{NULL} to omit. Defaults to
#'   \code{"I(julian_std^2)"}.
#' @param n_pred Integer. Number of prediction points. Defaults to \code{100L}.
#' @param nsamp Integer. Number of posterior samples.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{pred_det_df}}{Matrix/data frame with columns
#'     \code{pred_julian}, \code{seq_julian}, \code{mean},
#'     \code{criL}, \code{criU}.}
#'   \item{\code{seq_julian}}{Raw Julian day sequence.}
#'   \item{\code{pred_julian}}{Scaled Julian day sequence.}
#' }
#'
#' @export
predict_detection_julian <- function(out_bestmodel,
                                     julian_day,
                                     m_julian,
                                     sd_julian,
                                     linear_coef    = "julian_std",
                                     quadratic_coef = "I(julian_std^2)",
                                     n_pred         = 100L,
                                     nsamp) {

  stopifnot(
    "julian_day must be a matrix or data frame" =
      is.matrix(julian_day) || is.data.frame(julian_day),
    "nsamp must be a single positive integer" =
      length(nsamp) == 1L && nsamp > 0L
  )

  # Raw and scaled sequences
  seq_julian  <- seq(min(julian_day, na.rm = TRUE),
                     max(julian_day, na.rm = TRUE),
                     length.out = n_pred)
  pred_julian <- (seq_julian - m_julian) / sd_julian

  # Fill prediction array
  predDet <- array(NA_real_, dim = c(n_pred, nsamp))

  for (i in seq_len(nsamp)) {

    linear_effect <- out_bestmodel$alpha.comm.samples[i, linear_coef] *
      pred_julian

    quad_effect <- if (!is.null(quadratic_coef)) {
      out_bestmodel$alpha.comm.samples[i, quadratic_coef] * pred_julian^2
    } else {
      0
    }

    predDet[, i] <- plogis(
      out_bestmodel$alpha.comm.samples[i, "(Intercept)"] +
        linear_effect +
        quad_effect
    )
  }

  pred_det_df <- cbind(
    pred_julian = pred_julian,
    seq_julian  = seq_julian,
    mean        = apply(predDet, 1, mean),
    criL        = apply(predDet, 1, function(x) quantile(x, prob = 0.025)),
    criU        = apply(predDet, 1, function(x) quantile(x, prob = 0.975))
  )

  list(
    pred_det_df = pred_det_df,
    seq_julian  = seq_julian,
    pred_julian = pred_julian
  )
}


# =============================================================================
# 10. plot_detection_julian
# =============================================================================

#' Plot Community-Level Detection Probability Over Survey Dates
#'
#' Converts a scaled Julian day sequence to calendar dates and produces a
#' ribbon + line plot of community-level detection probability.
#'
#' @param pred_det_df Matrix or data frame as returned in
#'   \code{\link{predict_detection_julian}$pred_det_df}. Must have columns
#'   \code{seq_julian}, \code{mean}, \code{criL}, \code{criU}.
#' @param start_date Date or character (coercible to \code{Date}). Origin date
#'   used to convert Julian days to calendar dates.
#' @param end_date Date or character. End date for the x-axis scale.
#' @param color_ribbon Character. Ribbon fill colour. Defaults to
#'   \code{"#CCCCCC"}.
#' @param color_line Character. Line colour. Defaults to \code{"#1B3A5C"}.
#' @param output_path Character. If provided, exported as PNG.
#'   Defaults to \code{NULL}.
#' @param width,height Numeric. PNG dimensions in cm.
#' @param dpi Integer. Resolution.
#' @param ribbon_alpha Numeric. Transparency of the credible interval ribbon.
#'   Defaults to \code{0.55}.
#'
#' @return Invisibly returns the \code{ggplot} object.
#'
#' @export

plot_detection_julian <- function(pred_det_df,
                                  start_date,
                                  end_date,
                                  color_ribbon  = "#CCCCCC",
                                  color_line    = "#1B3A5C",
                                  ribbon_alpha  = 0.55,
                                  output_path   = NULL,
                                  width         = 28,
                                  height        = 24,
                                  dpi           = 300) {

  start_date <- as.Date(start_date)
  end_date   <- as.Date(end_date)

  month_sampling <- as.Date(pred_det_df[, "seq_julian"],
                            origin = start_date)

  p <- pred_det_df %>%
    tibble::as_tibble() %>%
    ggplot2::ggplot() +

    # # Subtle horizontal reference line at 0.5
    # ggplot2::geom_hline(
    #   yintercept = 0.5,
    #   linetype   = "dashed",
    #   colour     = "grey60",
    #   linewidth  = 0.4
    # ) +

    # 95% credible interval ribbon — semi-transparent so the line
    # remains visible even when the CI is very wide
    ggplot2::geom_ribbon(
      ggplot2::aes(x = month_sampling, ymin = criL, ymax = criU),
      fill  = color_ribbon,
      alpha = ribbon_alpha
    ) +

    # Mean prediction line — slightly thicker and on top of the ribbon
    ggplot2::geom_line(
      ggplot2::aes(x = month_sampling, y = mean),
      linewidth = 2.2,
      color     = color_line
    ) +

    # Axis labels
    ggplot2::labs(
      x        = "Survey Date (derived from Julian Day)",
      y        = "Detection probability"
    ) +

    # # Fix y-axis to the probability scale
    # ggplot2::scale_y_continuous(
    #   limits = c(0, 1),
    #   breaks = seq(0, 1, by = 0.2),
    #   expand = ggplot2::expansion(mult = c(0.01, 0.02))
    # ) +

    ggplot2::scale_x_date(
      date_labels       = "%b-%d",
      breaks            = seq(start_date, end_date, by = "2 weeks"),
      date_minor_breaks = "1 week",
      expand            = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +

    ggplot2::theme_classic() +
    ggplot2::theme(
      # Axis titles — slightly smaller than before for better balance
      axis.title.x  = ggplot2::element_text(size = 18, color = "#060508",
                                            margin = ggplot2::margin(t = 10)),
      axis.title.y  = ggplot2::element_text(size = 18, color = "#060508",
                                            margin = ggplot2::margin(r = 10)),
      # Axis tick labels
      axis.text.x   = ggplot2::element_text(color = "grey20", size = 13,
                                            angle = 45, vjust = 1, hjust = 1),
      axis.text.y   = ggplot2::element_text(color = "grey20", size = 14),
      # Thin axis lines
      axis.line     = ggplot2::element_line(colour = "grey70", linewidth = 0.4),
      axis.ticks    = ggplot2::element_line(colour = "grey70", linewidth = 0.4),
      # Minimal panel — only horizontal major gridlines for readability
      panel.grid.major.y = ggplot2::element_line(colour = "grey92",
                                                  linewidth = 0.3),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      # Caption styling
      plot.caption  = ggplot2::element_text(size = 11, colour = "grey50",
                                            hjust = 0,
                                            margin = ggplot2::margin(t = 8)),
      # Breathing room around the plot
      plot.margin   = ggplot2::margin(t = 16, r = 20, b = 10, l = 10)
    )

  if (!is.null(output_path)) {
    png(output_path, width = width, height = height,
        units = "cm", res = dpi, pointsize = 6, bg = "white")
    print(p)
    dev.off()
    message("Detection plot saved: ", basename(output_path))
  } else {
    print(p)
  }

  return(p)
}



# =============================================================================
# 11. prepare_occ_sequences
# =============================================================================

#' Prepare Raw and Scaled Prediction Sequences per Covariate
#'
#' For each continuous variable, generates an evenly-spaced sequence of
#' \code{n_pred} values across the observed range plus a mean-centred,
#' scaled version ready for model prediction.
#'
#' @param variables_to_scale Character vector of continuous variable names.
#'   Must match columns in \code{variables_df}.
#' @param variables_df Data frame containing the raw covariate data.
#' @param n_pred Integer. Number of prediction points. Defaults to \code{100L}.
#'
#' @return Named list with \code{seq_list} (raw) and
#'   \code{scaled_centered_list} (scaled), each a named list of numeric
#'   vectors of length \code{n_pred}.
#'
#' @export
prepare_occ_sequences <- function(variables_to_scale,
                                   variables_df,
                                   n_pred = 100L) {

  stopifnot(
    "variables_to_scale must be a non-empty character vector" =
      is.character(variables_to_scale) && length(variables_to_scale) > 0,
    "variables_df must be a data frame" =
      is.data.frame(variables_df),
    "all variables_to_scale must be columns in variables_df" =
      all(variables_to_scale %in% names(variables_df)),
    "n_pred must be a single positive integer" =
      length(n_pred) == 1L && n_pred > 0L
  )

  seq_list             <- vector("list", length(variables_to_scale))
  scaled_centered_list <- vector("list", length(variables_to_scale))
  names(seq_list)             <- variables_to_scale
  names(scaled_centered_list) <- variables_to_scale

  for (variable in variables_to_scale) {
    col      <- variables_df[[variable]]
    seq_vals  <- seq(min(col, na.rm = TRUE), max(col, na.rm = TRUE),
                     length.out = n_pred)
    pred_vals <- (seq_vals - mean(col, na.rm = TRUE)) / sd(col, na.rm = TRUE)
    seq_list[[variable]]             <- seq_vals
    scaled_centered_list[[variable]] <- pred_vals
  }

  list(seq_list = seq_list, scaled_centered_list = scaled_centered_list)
}


# =============================================================================
# 12. predict_occupancy
# =============================================================================

#' Fill Occupancy Prediction Array from MCMC Samples
#'
#' Initialises a \code{[n_pred × nsamp × n_combinations]} array and fills it
#' using community-level beta samples, iterating over variables, years, and
#' (optionally) type levels.
#'
#' @param out_bestmodel Fitted model with \code{beta.comm.samples}.
#' @param scaled_centered_list Named list of scaled sequences from
#'   \code{\link{prepare_occ_sequences}}.
#' @param variables_to_scale Character vector of continuous variable names.
#' @param year_levels Character vector of year/period levels as they appear
#'   in \code{beta.comm.samples} (use \code{"(Intercept)"} for single-year).
#' @param type_levels Character vector of type levels, or \code{NULL}.
#'   First element is treated as the reference (coefficient = 0).
#' @param nsamp Integer. Number of posterior samples.
#' @param n_pred Integer. Number of prediction points. Defaults to \code{100L}.
#'
#' @return Three-dimensional \code{NA_real_} array with named third dimension.
#'
#' @export
predict_occupancy <- function(out_bestmodel,
                               scaled_centered_list,
                               variables_to_scale,
                               year_levels,
                               type_levels = NULL,
                               nsamp,
                               n_pred      = 100L) {

  stopifnot(
    "out_bestmodel must have beta.comm.samples" =
      !is.null(out_bestmodel$beta.comm.samples),
    "nsamp must be a single positive integer" =
      length(nsamp) == 1L && nsamp > 0L
  )

  n_vars  <- length(variables_to_scale)
  n_years <- length(year_levels)

  # Build third-dimension names
  if (is.null(type_levels)) {
    n_combinations   <- n_vars * n_years
    dimnames_predOcc <- character(n_combinations)
    for (vi in seq_len(n_vars))
      for (yi in seq_len(n_years)) {
        idx <- (vi - 1L) * n_years + yi
        dimnames_predOcc[idx] <- paste0(variables_to_scale[vi], "_",
                                        year_levels[yi])
      }
  } else {
    n_types          <- length(type_levels)
    n_combinations   <- n_vars * n_years * n_types
    dimnames_predOcc <- character(n_combinations)
    for (vi in seq_len(n_vars))
      for (yi in seq_len(n_years))
        for (ti in seq_len(n_types)) {
          idx <- (vi - 1L) * n_years * n_types +
                 (yi - 1L) * n_types + ti
          dimnames_predOcc[idx] <- paste0(variables_to_scale[vi], "_",
                                          year_levels[yi], "_",
                                          type_levels[ti])
        }
  }

  predOcc <- array(NA_real_, dim = c(n_pred, nsamp, n_combinations))
  dimnames(predOcc) <- list(NULL, NULL, dimnames_predOcc)

  for (i in seq_len(nsamp)) {
    idx <- 1L
    for (var in variables_to_scale) {
      pred_var <- scaled_centered_list[[var]]
      for (survey_year in year_levels) {
        year_eff <- ifelse(survey_year != "(Intercept)",
                           out_bestmodel$beta.comm.samples[i, survey_year], 0)
        if (is.null(type_levels)) {
          predOcc[, i, idx] <- plogis(
            out_bestmodel$beta.comm.samples[i, "(Intercept)"] +
              year_eff +
              out_bestmodel$beta.comm.samples[i, var] * pred_var)
          idx <- idx + 1L
        } else {
          for (type in type_levels) {
            type_eff <- ifelse(type != type_levels[1],
                               out_bestmodel$beta.comm.samples[i, type], 0)
            predOcc[, i, idx] <- plogis(
              out_bestmodel$beta.comm.samples[i, "(Intercept)"] +
                year_eff + type_eff +
                out_bestmodel$beta.comm.samples[i, var] * pred_var)
            idx <- idx + 1L
          }
        }
      }
    }
  }

  predOcc
}


# =============================================================================
# 13. summarise_occ_predictions
# =============================================================================

#' Summarise Occupancy Prediction Array into a Tidy Data Frame
#'
#' Converts the 3-D prediction array produced by \code{\link{predict_occupancy}}
#' into a tidy \code{tibble} with posterior mean, SD, and 95\% CRI per
#' combination of variable, year, and (optionally) class.
#'
#' @param predOcc Three-dimensional array \code{[n_pred x nsamp x combinations]}
#'   as returned by \code{\link{predict_occupancy}}.
#' @param year_levels Character vector of year/period levels.
#' @param type_levels Character vector of type levels, or \code{NULL}.
#' @param seq_list Named list of raw sequences from
#'   \code{\link{prepare_occ_sequences}}.
#' @param variable_labels Optional named character vector of display labels
#'   for each variable. If \code{NULL}, variable names are used as-is.
#'
#' @return A tidy \code{tibble} with columns \code{variable},
#'   \code{expanded_variable}, \code{year}, \code{class} (if applicable),
#'   \code{mean}, \code{sd}, \code{ciL}, \code{ciU}, and \code{seq_data}.
#'
#' @export
summarise_occ_predictions <- function(predOcc,
                                       year_levels,
                                       type_levels     = NULL,
                                       seq_list        = NULL,
                                       variable_labels = NULL) {

  stopifnot(
    "predOcc must be a 3D array" = is.array(predOcc) && length(dim(predOcc)) == 3
  )

  combo_names <- dimnames(predOcc)[[3]]
  n_pred      <- dim(predOcc)[1]

  rows <- lapply(combo_names, function(combo) {

    # Parse combo name back into parts
    parts <- strsplit(combo, "_")[[1]]

    if (is.null(type_levels)) {
      # Format: variable_year  (year is last part)
      year_val <- parts[length(parts)]
      var_val  <- paste(parts[-length(parts)], collapse = "_")
      class_val <- year_val
    } else {
      # Format: variable_year_type  (last two parts are year and type)
      type_val  <- parts[length(parts)]
      year_val  <- parts[length(parts) - 1]
      var_val   <- paste(parts[-c(length(parts) - 1, length(parts))],
                         collapse = "_")
      class_val <- paste0(year_val, "_", type_val)
    }

    slice <- predOcc[, , combo]

    tibble::tibble(
      variable          = var_val,
      expanded_variable = if (!is.null(variable_labels) &&
                               var_val %in% names(variable_labels))
                            variable_labels[[var_val]] else var_val,
      year              = year_val,
      class             = class_val,
      mean              = list(apply(slice, 1, mean)),
      sd                = list(apply(slice, 1, sd)),
      ciL               = list(apply(slice, 1, quantile, probs = 0.025)),
      ciU               = list(apply(slice, 1, quantile, probs = 0.975)),
      seq_data          = list(if (!is.null(seq_list) && var_val %in% names(seq_list))
                                 seq_list[[var_val]] else seq_len(n_pred))
    )
  })

  dplyr::bind_rows(rows)
}


# =============================================================================
# 14. run_occ_predictions  (master wrapper)
# =============================================================================

#' Run the Full Occupancy Prediction Pipeline
#'
#' Calls \code{\link{prepare_occ_sequences}}, \code{\link{predict_occupancy}},
#' and \code{\link{summarise_occ_predictions}} in sequence.
#'
#' @param variables_to_scale Character vector of continuous variable names.
#' @param variables_df Data frame with raw covariate data.
#' @param out_bestmodel Fitted model object.
#' @param year_levels_model Year/period levels as in \code{beta.comm.samples}.
#' @param type_levels Character vector or \code{NULL}.
#' @param nsamp Integer. Number of posterior samples.
#' @param n_pred Integer. Prediction points per variable. Defaults to \code{100L}.
#' @param variable_labels Optional named character vector of plot labels.
#'
#' @return Named list: \code{seq_list}, \code{scaled_centered_list},
#'   \code{predOcc}, \code{pred_occ_df}.
#'
#' @export
run_occ_predictions <- function(variables_to_scale,
                                 variables_df,
                                 out_bestmodel,
                                 year_levels_model,
                                 type_levels     = NULL,
                                 nsamp,
                                 n_pred          = 100L,
                                 variable_labels = NULL) {

  message("Step 1/3: Preparing sequences...")
  seqs <- prepare_occ_sequences(variables_to_scale, variables_df, n_pred)

  message("Step 2/3: Running prediction loop over ", nsamp, " samples...")
  predOcc <- predict_occupancy(
    out_bestmodel, seqs$scaled_centered_list, variables_to_scale,
    year_levels_model, type_levels, nsamp, n_pred
  )

  message("Step 3/3: Summarising...")
  pred_occ_df <- summarise_occ_predictions(
    predOcc, year_levels_model, type_levels, seqs$seq_list, variable_labels
  )

  message("Done. Array: ", paste(dim(predOcc), collapse = " x "))
  list(seq_list             = seqs$seq_list,
       scaled_centered_list = seqs$scaled_centered_list,
       predOcc              = predOcc,
       pred_occ_df          = pred_occ_df)
}


# =============================================================================
# 15. plot_occ_predictions
# =============================================================================

#' Build Occupancy Prediction Plot List
#'
#' Produces a named list of \code{ggplot} objects, one per combination of
#' variable × year × (optionally) class.
#'
#' @param pred_occ_df Tidy data frame from \code{\link{summarise_occ_predictions}}.
#' @param type_levels Character vector or \code{NULL}.
#' @param color_ribbon,color_line,color_text,color_strip_bg Characters.
#'   Colours for ribbon, line, text, and strip background.
#' @param strip_class_prefix Character to strip from class names in titles.
#' @param ribbon_alpha Numeric. Transparency of the credible interval ribbon.
#'   Defaults to \code{0.55}.
#'
#' @return Named list of \code{ggplot} objects.
#'
#' @export
plot_occ_predictions <- function(pred_occ_df,
                                  type_levels        = NULL,
                                  color_ribbon       = "#CCCCCC",
                                  color_line         = "#1B3A5C",
                                  color_text         = "#1B3A5C",
                                  color_strip_bg     = "#E8EEF4",
                                  ribbon_alpha       = 0.55,
                                  strip_class_prefix = "") {

  unique_years       <- unique(pred_occ_df$year)
  show_year_in_title <- length(unique_years) > 1

  grouping_vars <- if (is.null(type_levels)) {
    c("expanded_variable", "year")
  } else {
    c("expanded_variable", "year", "class")
  }

  make_title <- if (is.null(type_levels)) {
    function(d) if (show_year_in_title) paste("Survey Year:", d$year[1]) else ""
  } else {
    function(d) {
      cl <- gsub(strip_class_prefix, "", d$class[1])
      if (show_year_in_title) paste("Land use:", cl, "- Year:", d$year[1])
      else                     paste("Land use:", cl)
    }
  }

  splits <- pred_occ_df %>%
    dplyr::group_split(!!!rlang::syms(grouping_vars))

  plot_names <- sapply(splits, function(d) {
    if (is.null(type_levels)) paste0(d$variable[1], "_", d$year[1])
    else                       paste0(d$variable[1], "_", d$year[1], "_", d$class[1])
  })

  plot_list <- lapply(splits, function(data_subset) {
    data_subset <- data_subset %>%
      tidyr::unnest(cols = c(mean, sd, ciL, ciU, seq_data))

    cat("Plotting:", data_subset$variable[1], "| Year:", data_subset$year[1])
    if (!is.null(type_levels)) cat(" | Class:", data_subset$class[1])
    cat(" |", nrow(data_subset), "rows\n")

    title_text <- make_title(data_subset)

    ggplot2::ggplot(data_subset,
                    ggplot2::aes(x = as.numeric(seq_data), y = mean)) +

      # # Subtle reference line at 0.5
      # ggplot2::geom_hline(
      #   yintercept = 0.5,
      #   linetype   = "dashed",
      #   colour     = "grey60",
      #   linewidth  = 0.4
      # ) +

      # 95% credible interval ribbon — semi-transparent
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = ciL, ymax = ciU),
        fill  = color_ribbon,
        alpha = ribbon_alpha
      ) +

      # Mean prediction line on top
      ggplot2::geom_line(
        color     = color_line,
        linewidth = 2.2
      ) +

      ggplot2::labs(
        x       = data_subset$expanded_variable[1],
        y       = "Occupancy probability (Community-level)",
        title   = title_text,
        caption = "Shaded area: 95% credible interval"
      ) +

      # # Fix y-axis to the full probability scale
      # ggplot2::scale_y_continuous(
      #   limits = c(0, 1),
      #   breaks = seq(0, 1, by = 0.2),
      #   expand = ggplot2::expansion(mult = c(0.01, 0.02))
      # ) +

      # ggplot2::scale_x_continuous(
      #   expand = ggplot2::expansion(mult = c(0.01, 0.01))
      # ) +

      ggplot2::theme_classic() +
      ggplot2::theme(

        # Axis titles
        axis.title.x  = ggplot2::element_text(size = 14, color = color_text,
                                              margin = ggplot2::margin(t = 10)),
        axis.title.y  = ggplot2::element_text(size = 14, color = color_text,
                                              margin = ggplot2::margin(r = 10)),

        # Axis tick labels
        axis.text     = ggplot2::element_text(size = 12, color = "grey20"),

        # Thin axis lines and ticks
        axis.line     = ggplot2::element_line(colour = "grey70", linewidth = 0.4),
        axis.ticks    = ggplot2::element_line(colour = "grey70", linewidth = 0.4),

        # Only horizontal gridlines for reading probability values
        panel.grid.major.y = ggplot2::element_line(colour = "grey92",
                                                    linewidth = 0.3),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),

        # Plot title (shown only when multiple years or classes)
        plot.title    = ggplot2::element_text(size = 13, color = color_text,
                                              face = "bold",
                                              margin = ggplot2::margin(b = 6)),

        # Caption
        plot.caption  = ggplot2::element_text(size = 10, colour = "grey50",
                                              hjust = 0,
                                              margin = ggplot2::margin(t = 8)),

        # Strip labels (for faceted use if ever needed)
        strip.text.x     = ggplot2::element_text(size = 12, color = color_text),
        strip.background = ggplot2::element_rect(fill = color_strip_bg,
                                                 colour = NA),

        # Breathing room
        plot.margin   = ggplot2::margin(t = 16, r = 20, b = 10, l = 10),
        legend.position = "none"
      )
  })

  names(plot_list) <- plot_names
  plot_list
}



# =============================================================================
# 16. export_occ_plots
# =============================================================================

#' Export Occupancy Prediction Plots to PNG
#'
#' Saves individual PNGs for every plot in \code{plot_list} and optionally
#' assembles a multiplot for a selected subset of variables.
#'
#' @param plot_list Named list of \code{ggplot} objects from
#'   \code{\link{plot_occ_predictions}}.
#' @param output_dir Character. Directory for output files (trailing \code{/}
#'   added if missing).
#' @param project_name Character. Appended to all filenames.
#' @param type_levels Character vector or \code{NULL} (affects filename).
#' @param selected_vars Character vector of variable names for the multiplot,
#'   or \code{NULL} to include all.
#' @param width_single,height_single Numeric. Individual plot dimensions in cm.
#' @param width_multi,height_multi Numeric. Multiplot dimensions in cm.
#' @param dpi Integer. Resolution.
#'
#' @return Invisibly returns the multiplot or \code{NULL}.
#'
#' @export
export_occ_plots <- function(plot_list,
                              output_dir,
                              project_name,
                              type_levels    = NULL,
                              selected_vars  = NULL,
                              width_single   = 24,
                              height_single  = 18,
                              width_multi    = 30,
                              height_multi   = 20,
                              dpi            = 600) {

  if (!endsWith(output_dir, "/")) output_dir <- paste0(output_dir, "/")

  # Individual plots
  for (p in plot_list) {
    variable_name <- p$data$variable[1]
    filename <- paste0(output_dir, "community_Occupancy_",
                       variable_name, "_", project_name, ".png")
    png(filename, width = width_single, height = height_single,
        units = "cm", res = dpi, bg = "white")
    print(p)
    dev.off()
    cat("Exported:", basename(filename), "\n")
  }

  # Multiplot
  sig_plots <- if (!is.null(selected_vars)) {
    Filter(function(p) p$data$variable[1] %in% selected_vars, plot_list)
  } else {
    plot_list
  }

  if (length(sig_plots) < 2) {
    cat("Fewer than 2 variables for multiplot — skipping.\n")
    return(invisible(NULL))
  }

  n_plots <- length(sig_plots)
  ncol_mp <- if (n_plots <= 4) 2 else if (n_plots <= 9) 3 else 4
  nrow_mp <- ceiling(n_plots / ncol_mp)

  x_labels <- sapply(sig_plots, function(p) p$labels$x)
  same_x   <- length(unique(x_labels)) == 1

  sig_plots <- lapply(sig_plots, function(p) {
    if (same_x) p + ggplot2::labs(y = NULL, x = NULL)
    else        p + ggplot2::labs(y = NULL)
  })

  unique_years <- unique(sapply(plot_list, function(p) p$data$year[1]))

  mf <- paste0(
    output_dir,
    "community_Occupancy_significant_multiplot_",
    if (!is.null(type_levels)) "classes_" else "",
    if (length(unique_years) > 1) paste0(paste(unique_years, collapse = "_"), "_") else "",
    project_name, ".png"
  )

  multiplot <- ggpubr::ggarrange(
    plotlist   = sig_plots,
    ncol       = ncol_mp, nrow = nrow_mp,
    labels     = "AUTO",
    font.label = list(size = 14, face = "bold")
  ) %>%
    ggpubr::annotate_figure(
      left   = ggpubr::text_grob("Occupancy probability (Community-level)",
                                 rot = 90, size = 14, face = "bold"),
      bottom = if (same_x) ggpubr::text_grob(x_labels[1], size = 14,
                                              face = "bold") else NULL
    )

  ggplot2::ggsave(mf, multiplot, width = width_multi, height = height_multi,
                  units = "cm", dpi = dpi, bg = "white")
  cat("Exported multiplot:", basename(mf), "\n")
  invisible(multiplot)
}


# =============================================================================
# 17. extract_psi_stats
# =============================================================================

#' Extract Site-Level Occupancy Probability Statistics
#'
#' Computes mean, SD, and 95\% credible interval bounds from the
#' three-dimensional \code{psi.samples} array (posterior samples ×
#' species × sites).
#'
#' @param out_bestmodel Fitted model with \code{psi.samples},
#'   \code{sp.names}, and associated site names in
#'   \code{data.list$occ.covs$sites}.
#' @param site_names Character vector of site names to use as column labels.
#'
#' @return Named list with four matrices (\code{mean}, \code{sd},
#'   \code{cri025}, \code{cri975}), each with species as rows and sites
#'   as columns.
#'
#' @export
extract_psi_stats <- function(out_bestmodel, site_names) {

  stopifnot(
    "out_bestmodel must have psi.samples" =
      !is.null(out_bestmodel$psi.samples),
    "site_names must match the site dimension of psi.samples" =
      length(site_names) == dim(out_bestmodel$psi.samples)[3]
  )

  sp_names <- out_bestmodel$sp.names

  add_dimnames <- function(mat) {
    rownames(mat) <- sp_names
    colnames(mat) <- site_names
    mat
  }

  list(
    mean   = add_dimnames(apply(out_bestmodel$psi.samples, c(2, 3), mean)),
    sd     = add_dimnames(apply(out_bestmodel$psi.samples, c(2, 3), sd)),
    cri025 = add_dimnames(apply(out_bestmodel$psi.samples, c(2, 3),
                                function(x) quantile(x, prob = 0.025))),
    cri975 = add_dimnames(apply(out_bestmodel$psi.samples, c(2, 3),
                                function(x) quantile(x, prob = 0.975)))
  )
}


# =============================================================================
# 18. export_psi_tables
# =============================================================================

#' Export Site-Level Occupancy Probability Tables to CSV
#'
#' Writes one CSV per statistic (mean, SD, CRI 0.025, CRI 0.975) from the
#' list returned by \code{\link{extract_psi_stats}}. Tables are transposed to
#' sites × species format with a leading \code{sites} column.
#'
#' @param psi_stats Named list as returned by \code{\link{extract_psi_stats}}.
#' @param output_dir Character. Output directory (trailing \code{/} added if
#'   missing).
#' @param project_name Character. Appended to all filenames.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect.
#'
#' @export
export_psi_tables <- function(psi_stats, output_dir, project_name) {

  stopifnot(
    "psi_stats must be a named list" =
      is.list(psi_stats) && !is.null(names(psi_stats))
  )

  if (!endsWith(output_dir, "/")) output_dir <- paste0(output_dir, "/")

  stat_names <- list(
    mean   = "psi_mean",
    sd     = "psi_sd",
    cri025 = "psi_cri025",
    cri975 = "psi_cri975"
  )

  for (stat in names(stat_names)) {

    filepath <- paste0(output_dir, stat_names[[stat]], "_",
                       project_name, ".csv")

    psi_stats[[stat]] %>%
      t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "sites") %>%
      readr::write_csv(filepath)

    message("Saved: ", basename(filepath))
  }

  invisible(NULL)
}
