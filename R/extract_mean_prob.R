#' Extract Species-Level Occupancy or Detection Probabilities
#'
#' Back-transforms MCMC samples to the probability scale (logit → probability)
#' and summarises species-level estimates across one or more classes.
#'
#' @author Thiago Bicudo
#'
#' @param samples Matrix of MCMC samples (\code{[nsamp × parameters]}).
#'   For occupancy use \code{out_bestmodel$beta.samples}; for detection use
#'   \code{out_bestmodel$alpha.samples}.
#'
#' @param class_list Named list mapping class labels to coefficient names,
#'   or \code{NULL} for intercept-only models.
#'
#'   Example (intercept only):
#'   \code{list(survey2025 = NULL)}
#'
#'   Example (categorical effect with reference level):
#'   \code{list(survey2022 = NULL, survey2023 = "year2023")}
#'
#' @param spp_reference Data frame with at least a \code{species} column,
#'   used to validate and join species names.
#'
#' @return A \code{tibble} with columns:
#'   \code{species}, \code{param_spp}, \code{class},
#'   \code{mean}, \code{sd}, \code{ciL}, \code{ciU}.
#'
#' @importFrom dplyr bind_rows mutate left_join
#' @importFrom tibble as_tibble
#'
#' @export
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
    c(
      mean = mean(x),
      sd   = sd(x),
      ciL  = quantile(x, 0.025),
      ciU  = quantile(x, 0.975)
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
