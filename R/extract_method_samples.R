#' Extract alpha.samples columns for a specific detection method
#'
#' Slices the correct columns from \code{alpha.samples} for a given detection
#' method in an integrated multi-species occupancy model, based on the number
#' of species and parameters per method. Column blocks are assumed to follow
#' the \code{spOccupancy} convention: all parameters for method 1 precede all
#' parameters for method 2, with each block structured as
#' \code{[param1_spp1, ..., param1_sppN, param2_spp1, ..., paramK_sppN]}.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @param alpha_samples A matrix or \code{mcmc} object of MCMC samples with
#'   named columns following the \code{spOccupancy} convention
#'   \code{"parameter-species"} (e.g., \code{out_bestmodel$alpha.samples}).
#' @param method_spec A named list describing each detection method, in the
#'   same order as they appear in \code{alpha.samples}. Each element must be
#'   a list with two entries:
#'   \describe{
#'     \item{\code{n_spp}}{Integer. Number of species for this method
#'       (e.g., \code{nrow(data.list$y$audio)}).}
#'     \item{\code{n_params}}{Integer. Number of parameters per species,
#'       including the intercept
#'       (e.g., \code{ncol(out_bestmodel$X.p[[1]])}).}
#'   }
#' @param target Character. Name of the method to extract, matching a name
#'   in \code{method_spec}.
#'
#' @return A matrix of dimensions \code{[nsamp x (n_spp * n_params)]} with
#'   the columns of \code{alpha_samples} corresponding to the target method.
#'   Column names are preserved from the input.
#'
#' @details
#' In integrated models fitted with \code{spOccupancy}, \code{alpha.samples}
#' contains species-specific detection coefficients for all methods
#' concatenated into a single matrix. Each method occupies a contiguous block
#' of columns, ordered as: all parameters for method 1, then all parameters
#' for method 2, and so on.
#'
#' Because column names are frequently duplicated across methods (e.g., two
#' \code{"(Intercept)"} columns — one per method), selecting columns by name
#' would silently return incorrect results. This function instead computes
#' block boundaries from cumulative column offsets derived from
#' \code{method_spec}, ensuring correct and unambiguous slicing by position.
#'
#' The \code{method_spec} list must include ALL methods in the model and in
#' the correct order. Changing the order will produce incorrect slices without
#' raising an error. Use the following to verify the order and content before
#' defining \code{method_spec}:
#'
#' \preformatted{
#' # Inspect parameters per method
#' sapply(out_bestmodel$X.p, colnames)
#'
#' # Build method_spec from model output
#' method_spec <- list(
#'   audio = list(
#'     n_spp    = nrow(data.list$y$audio),
#'     n_params = ncol(out_bestmodel$X.p[[1]])
#'   ),
#'   camtrap = list(
#'     n_spp    = nrow(data.list$y$camtrap),
#'     n_params = ncol(out_bestmodel$X.p[[2]])
#'   )
#' )
#'
#' # Sanity check: block sizes must sum to ncol(alpha.samples)
#' stopifnot(
#'   ncol(as.matrix(out_bestmodel$alpha.samples)) ==
#'     nrow(data.list$y$audio)   * ncol(out_bestmodel$X.p[[1]]) +
#'     nrow(data.list$y$camtrap) * ncol(out_bestmodel$X.p[[2]])
#' )
#'
#' # Extract per method
#' samples_aud <- extract_method_samples(
#'   alpha_samples = out_bestmodel$alpha.samples,
#'   method_spec   = method_spec,
#'   target        = "audio"
#' )
#'
#' samples_cam <- extract_method_samples(
#'   alpha_samples = out_bestmodel$alpha.samples,
#'   method_spec   = method_spec,
#'   target        = "camtrap"
#' )
#' }
#'
#' @export
extract_method_samples <- function(alpha_samples, method_spec, target) {

  # ── Input validation ───────────────────────────────────────────────────────
  if (inherits(alpha_samples, "mcmc")) {
    alpha_samples <- as.matrix(alpha_samples)
  }

  stopifnot(
    "alpha_samples must be a matrix or mcmc object" =
      is.matrix(alpha_samples),
    "alpha_samples must have named columns"         =
      !is.null(colnames(alpha_samples)),
    "method_spec must be a named list"              =
      is.list(method_spec) && !is.null(names(method_spec)),
    "each element of method_spec must contain n_spp and n_params" =
      all(sapply(method_spec, function(m) all(c("n_spp", "n_params") %in% names(m)))),
    "target must be a single character string"      =
      is.character(target) && length(target) == 1,
    "target must match a name in method_spec"       =
      target %in% names(method_spec)
  )

  # ── Compute cumulative column offsets ──────────────────────────────────────
  block_sizes <- sapply(method_spec, function(m) m$n_spp * m$n_params)
  cum_end     <- cumsum(block_sizes)
  cum_start   <- c(1, cum_end[-length(cum_end)] + 1)
  names(cum_start) <- names(method_spec)
  names(cum_end)   <- names(method_spec)

  # ── Validate total against actual dimensions ───────────────────────────────
  total_cols <- sum(block_sizes)

  if (total_cols != ncol(alpha_samples)) {
    stop(
      "Total columns derived from method_spec (", total_cols, ") does not ",
      "match ncol(alpha_samples) (", ncol(alpha_samples), ").\n",
      "Check that n_spp and n_params are correct for each method and that ",
      "all methods are listed in method_spec."
    )
  }

  # ── Slice target block ─────────────────────────────────────────────────────
  idx <- cum_start[target]:cum_end[target]
  alpha_samples[, idx, drop = FALSE]
}
