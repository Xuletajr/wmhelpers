#' Compute mean detection or occupancy probability per species from MCMC samples
#'
#' Back-transforms MCMC samples to the probability scale via \code{plogis()}
#' and returns posterior summaries of detection or occupancy probability for
#' each species, optionally adding species-specific covariate coefficients to
#' the intercept. Designed for outputs of multi-species occupancy models fitted
#' with \code{spOccupancy}, including single-method and integrated models.
#'
#' @author Thiago Bicudo
#'
#' @param samples A matrix or \code{mcmc} object of MCMC samples with named
#'   columns following the pattern \code{"parameter-species"} (e.g.,
#'   \code{"(Intercept)-Alethe castanea"}, \code{"numD-Alethe castanea"}).
#'   For occupancy use \code{out_bestmodel$beta.samples}; for detection use
#'   \code{out_bestmodel$alpha.samples}. For integrated models, pass the
#'   relevant column slice for each data source separately
#'   (e.g., \code{as.matrix(alpha.samples)[, 1:267]}).
#'
#' @param class_list A named list where each element defines one category.
#'   Use \code{NULL} for the reference category (intercept only). For other
#'   categories, provide the coefficient prefix as a character string matching
#'   the column naming pattern in \code{samples}.
#'
#'   Example (intercept only):
#'   \code{list(survey2025 = NULL)}
#'
#'   Example (categorical effect with reference level):
#'   \code{list(survey2022 = NULL, survey2023 = "year2023")}
#'
#' @param spp_reference A data frame with a column named \code{species},
#'   used to validate species names against those derived from column names.
#'
#' @return A \code{tibble} with one row per species per class, containing
#'   columns \code{species}, \code{param_spp}, \code{class}, \code{mean},
#'   \code{sd}, \code{ciL}, and \code{ciU}.
#'
#' @details
#' Column names in \code{samples} must follow the \code{spOccupancy} convention
#' \code{"parameter-species"}. Intercept columns are identified by the pattern
#' \code{"(Intercept)"}. Species names are extracted from these columns and used
#' to construct the full column names for any additional coefficient
#' (e.g., \code{"year2023-Alethe castanea"}).
#'
#' For integrated models (e.g., audio + camera trap), call this function
#' separately for each data source by slicing the relevant columns from
#' \code{alpha.samples} before passing to this function.
#'
#' A warning is issued if the species derived from column names differ in
#' composition or order from \code{spp_reference}. Execution proceeds using
#' the order from column names, which is the order used internally for all
#' computations.
#'
#' @importFrom dplyr bind_rows as_tibble select
#'
#' @export
extract_mean_prob <- function(samples, class_list, spp_reference) {

  # ── Input validation ───────────────────────────────────────────────────────
  if (inherits(samples, "mcmc")) {
    samples <- as.matrix(samples)
  }

  stopifnot(
    "samples must be a matrix or mcmc object"      = is.matrix(samples),
    "samples must have named columns"              = !is.null(colnames(samples)),
    "class_list must be a named list"              = is.list(class_list) &&
      !is.null(names(class_list)),
    "spp_reference must be a data frame"           = is.data.frame(spp_reference),
    "spp_reference must contain a 'species' column" = "species" %in%
      names(spp_reference)
  )

  # ── Extract intercept columns and derive species names ────────────────────
  intercept_cols <- grep("\\(Intercept\\)", colnames(samples), value = TRUE)

  if (length(intercept_cols) == 0) {
    stop("No '(Intercept)' columns found in samples.")
  }

  spp_names <- sub("^[^-]+-", "", intercept_cols)

  if (!identical(spp_names, spp_reference$species)) {
    warning(
      "Species derived from intercept columns do not match spp_reference. ",
      "Proceeding with species names extracted from column names."
    )
  }

  # ── Posterior summary helper ───────────────────────────────────────────────
  summarise_post <- function(x) {
    q <- quantile(x, c(0.025, 0.975), names = FALSE)
    c(mean = mean(x), sd = sd(x), ciL = q[1], ciU = q[2])
  }

  # ── Compute per class ──────────────────────────────────────────────────────
  results <- lapply(names(class_list), function(cls) {

    coef_prefix <- class_list[[cls]]

    # Reference category: intercept only
    linpred <- samples[, intercept_cols, drop = FALSE]  # [nsamp x n_spp]

    # Non-reference category: add species-specific coefficient
    if (!is.null(coef_prefix)) {

      coef_cols <- paste0(coef_prefix, "-", spp_names)
      missing   <- setdiff(coef_cols, colnames(samples))

      if (length(missing) > 0) {
        stop(
          "The following coefficient columns were not found in samples ",
          "for prefix '", coef_prefix, "':\n",
          paste(missing, collapse = "\n"),
          "\nCheck that class_list prefixes match the column naming pattern ",
          "in samples (expected: '<prefix>-<species>')."
        )
      }

      coef_mat <- samples[, coef_cols, drop = FALSE]  # [nsamp x n_spp]
      linpred  <- linpred + coef_mat
    }

    # Probability scale and posterior summary
    prob  <- plogis(linpred)
    stats <- apply(prob, 2, summarise_post)

    df           <- as.data.frame(t(stats))
    df$param_spp <- if (is.null(coef_prefix)) intercept_cols else coef_cols
    df$species   <- spp_names
    df$class     <- cls
    df
  })

  # ── Bind and format output ─────────────────────────────────────────────────
  out <- dplyr::bind_rows(results)

  dplyr::as_tibble(out) |>
    dplyr::select(species, param_spp, class, mean, sd, ciL, ciU)
}

extract_mean_prob <- function(samples, class_list, spp_reference) {

  # ── Validation ──────────────────────────────────────────────────────────────
  stopifnot(
    is.matrix(samples),
    !is.null(colnames(samples)),
    is.list(class_list),
    !is.null(names(class_list)),
    is.data.frame(spp_reference),
    "species" %in% names(spp_reference)
  )

  # ── Identify intercept columns (species-level intercepts) ───────────────────
  intercept_cols <- grep("\\(Intercept\\)", colnames(samples), value = TRUE)

  if (length(intercept_cols) == 0) {
    stop("No '(Intercept)' columns found in samples.")
  }

  # Extract species names from parameter names
  spp_names <- gsub(".*\\]", "", intercept_cols)  # fallback safe
  # Better approach:
  spp_names <- gsub(".*species\\[(\\d+)\\]", "\\1", intercept_cols)

  # ── Helper: summarise posterior ─────────────────────────────────────────────
  summarise_post <- function(x) {
    q <- quantile(x, c(0.025, 0.975), names = FALSE)

    c(
      mean = mean(x),
      sd   = sd(x),
      ciL  = q[1],
      ciU  = q[2]
    )
  }

  # ── Loop over classes ───────────────────────────────────────────────────────
  results <- lapply(names(class_list), function(cls) {

    coef_name <- class_list[[cls]]

    # Base: intercept
    linpred <- samples[, intercept_cols, drop = FALSE]

    # Add class effect if exists
    if (!is.null(coef_name)) {

      if (!coef_name %in% colnames(samples)) {
        stop(paste("Coefficient", coef_name, "not found in samples"))
      }

      linpred <- linpred + samples[, coef_name]
    }

    # Back-transform
    prob <- plogis(linpred)

    # Summarise per species (column-wise)
    stats <- apply(prob, 2, summarise_post)

    df <- as.data.frame(t(stats))
    df$param_spp <- intercept_cols
    df$class <- cls

    return(df)
  })

  # ── Combine ─────────────────────────────────────────────────────────────────
  out <- dplyr::bind_rows(results)

  # ── Attach species names ────────────────────────────────────────────────────
  out$species <- spp_reference$species

  # ── Final formatting ────────────────────────────────────────────────────────
  out <- dplyr::as_tibble(out) %>%
    dplyr::select(species, param_spp, class, mean, sd, ciL, ciU)

  return(out)
}
