# =============================================================================
# Jaccard similarity computation and visualization functions
# File: R/jaccard_index.R
# =============================================================================


#' Compute pairwise Jaccard similarity from posterior Z samples
#'
#' Computes pairwise species composition similarity (Jaccard index) across
#' MCMC samples from the latent occupancy matrix (Z). Supports site-level
#' comparisons or group-level aggregation (e.g., by region), where a species
#' is considered present in a group if detected in any site within that group
#' for a given MCMC sample.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta") and Thiago Bicudo
#'
#' @param z_samples A 3-dimensional array \code{[MCMC x species x sites]},
#'   typically \code{out_bestmodel$z.samples}.
#' @param site_names Character vector of site names aligned with the third
#'   dimension of \code{z_samples}.
#' @param group_var Optional named character vector mapping site names to group
#'   labels (e.g., region). Names must match \code{site_names}. If \code{NULL}
#'   (default), site-level Jaccard is computed. If provided, presence is
#'   aggregated across sites within each group per MCMC sample — a species is
#'   considered present in the group if detected in at least one site.
#' @param species_idx Optional integer vector to subset species before computing
#'   similarity. Use \code{get_species_idx()} to derive from a traits table.
#'   Defaults to \code{NULL}.
#'
#' @return A named list with three elements: \code{samples} (a matrix
#'   \code{[MCMC x n_pairs]} with per-sample Jaccard values), \code{summary}
#'   (a tibble with posterior summaries — mean, median, sd, cil, ciu,
#'   ci_width, cv — per pair), and \code{matrix_mean} (a square matrix of
#'   mean Jaccard values).
#'
#' @details
#' Jaccard similarity: \eqn{J = |A \cap B| / |A \cup B|}. Values near 1
#' indicate high compositional similarity; values near 0 indicate low overlap.
#' To convert to dissimilarity (beta diversity), use \code{dissimilarity = TRUE}
#' in \code{summarise_jaccard_groups()}.
#'
#' For group-level computation, a species is considered present in a group if
#' \code{z == 1} in at least one site within that group for a given MCMC
#' sample (\code{max} across sites). This is a union-based presence definition.
#'
#' Computation time scales with \code{n_mcmc x n_pairs}. For site-level
#' Jaccard with many sites (e.g., 150 sites = 22,500 pairs), expect several
#' minutes.
#'
#' @seealso \code{\link{summarise_jaccard_groups}}, \code{\link{plot_jaccard_heatmap}},
#'   \code{\link{plot_jaccard_groups}}
#'
#' @export
compute_jaccard <- function(z_samples,
                            site_names,
                            group_var   = NULL,
                            species_idx = NULL) {

  stopifnot(
    "z_samples must be a 3D array"         = is.array(z_samples) &&
      length(dim(z_samples)) == 3,
    "site_names must be a character vector" = is.character(site_names),
    "site_names length must match dim 3"    = length(site_names) ==
      dim(z_samples)[3]
  )

  if (!is.null(species_idx)) {
    z_samples <- z_samples[, species_idx, , drop = FALSE]
    message(length(species_idx), " species used for Jaccard computation.")
  }

  n_mcmc    <- dim(z_samples)[1]
  n_species <- dim(z_samples)[2]

  # ── Define comparison units (sites or groups) ─────────────────────────────
  if (!is.null(group_var)) {

    unit_names <- unique(group_var)
    unit_idx   <- lapply(unit_names, function(g) which(group_var == g))
    names(unit_idx) <- unit_names

    message("Group-level Jaccard: ", length(unit_names), " groups, ",
            length(unit_names)^2, " pairs.")

  } else {

    unit_names <- site_names
    unit_idx   <- lapply(seq_along(site_names), function(i) i)
    names(unit_idx) <- site_names

    message("Site-level Jaccard: ", length(site_names), " sites, ",
            length(site_names)^2, " pairs.")
  }

  n_units <- length(unit_names)

  # ── Pre-compute presence array [MCMC x species x unit] ────────────────────
  message("Pre-computing presence matrices...")

  presence_array <- array(NA_real_,
                          dim      = c(n_mcmc, n_species, n_units),
                          dimnames = list(NULL, NULL, unit_names))

  for (k in seq_len(n_mcmc)) {
    for (u in seq_along(unit_names)) {
      idx <- unit_idx[[u]]
      presence_array[k, , u] <- if (length(idx) == 1) {
        z_samples[k, , idx]
      } else {
        apply(z_samples[k, , idx], 1, max)
      }
    }
  }

  # ── Compute Jaccard for all pairs ─────────────────────────────────────────
  pairs      <- expand.grid(a = unit_names, b = unit_names,
                            stringsAsFactors = FALSE)
  pair_names <- paste0(pairs$a, "__", pairs$b)
  n_pairs    <- nrow(pairs)

  message("Computing Jaccard for ", n_mcmc, " MCMC samples x ",
          n_pairs, " pairs...")

  jacc_array           <- matrix(NA_real_, nrow = n_mcmc, ncol = n_pairs)
  colnames(jacc_array) <- pair_names

  for (k in seq_len(n_mcmc)) {
    for (p in seq_len(n_pairs)) {
      pa <- presence_array[k, , pairs$a[p]]
      pb <- presence_array[k, , pairs$b[p]]
      intersection     <- sum(pa * pb)
      union            <- sum(pa) + sum(pb) - intersection
      jacc_array[k, p] <- if (union > 0) intersection / union else NA_real_
    }
    if (k %% 500 == 0) message("  MCMC sample ", k, " of ", n_mcmc)
  }

  # ── Summarise posterior ────────────────────────────────────────────────────
  summary_df <- tibble::tibble(
    pair   = pair_names,
    unit_a = pairs$a,
    unit_b = pairs$b,
    mean   = apply(jacc_array, 2, mean,              na.rm = TRUE),
    median = apply(jacc_array, 2, stats::median,     na.rm = TRUE),
    sd     = apply(jacc_array, 2, stats::sd,         na.rm = TRUE),
    cil    = apply(jacc_array, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
    ciu    = apply(jacc_array, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  ) %>%
    dplyr::mutate(
      ci_width = ciu - cil,
      cv       = sd / mean
    )

  # ── Build square mean matrix ───────────────────────────────────────────────
  matrix_mean <- summary_df %>%
    dplyr::select(unit_a, unit_b, mean) %>%
    tidyr::pivot_wider(names_from = unit_b, values_from = mean) %>%
    tibble::column_to_rownames("unit_a") %>%
    as.matrix()

  matrix_mean <- matrix_mean[colnames(matrix_mean), ]

  list(
    samples     = jacc_array,
    summary     = summary_df,
    matrix_mean = matrix_mean
  )
}


#' Heatmap of pairwise Jaccard similarity
#'
#' Generates a tile heatmap of any posterior summary statistic from
#' \code{compute_jaccard()}, with flexible diagonal treatment and automatic
#' color scale adjustment by statistic type.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta") and Thiago Bicudo
#'
#' @param jaccard_result Output from \code{compute_jaccard()}.
#' @param title Character. Plot title. Defaults to
#'   \code{"Pairwise Jaccard Similarity"}.
#' @param stat Character. Posterior summary statistic to display. One of
#'   \code{"mean"}, \code{"median"}, \code{"sd"}, \code{"cil"}, \code{"ciu"},
#'   \code{"ci_width"}, or \code{"cv"}. Defaults to \code{"mean"}.
#' @param diagonal Character. Color string for diagonal tiles (e.g.,
#'   \code{"black"}, \code{"white"}, \code{"#1B3644"}), or \code{"value"} to
#'   apply the same color scale as off-diagonal tiles. Defaults to
#'   \code{"#1B3644"}.
#' @param show_values Logical. If \code{TRUE}, displays statistic values inside
#'   tiles. Defaults to \code{TRUE}.
#' @param digits Integer. Decimal places for tile labels. Defaults to \code{2}.
#' @param low_color Character. Hex color for the low end of the fill scale.
#'   Defaults to \code{"#BFFFE6"}.
#' @param high_color Character. Hex color for the high end of the fill scale.
#'   Defaults to \code{"#1B3644"}.
#' @param text_color Character. Hex color for axis text. Defaults to
#'   \code{"#1B3644"}.
#'
#' @return A \code{ggplot2} object.
#'
#' @details
#' Color scale limits are automatically adjusted by statistic type: similarity
#' statistics (\code{mean}, \code{median}, \code{cil}, \code{ciu}) use
#' \code{[0, 1]}; uncertainty statistics (\code{sd}, \code{ci_width},
#' \code{cv}) use \code{[0, observed max]}.
#'
#' @seealso \code{\link{compute_jaccard}}
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient labs theme_bw theme element_text element_blank geom_text
#'
#' @export
plot_jaccard_heatmap <- function(jaccard_result,
                                 title       = "Pairwise Jaccard Similarity",
                                 stat        = "mean",
                                 diagonal    = "#1B3644",
                                 show_values = TRUE,
                                 digits      = 2,
                                 low_color   = "#BFFFE6",
                                 high_color  = "#1B3644",
                                 text_color  = "#1B3644") {

  valid_stats <- c("mean", "median", "sd", "cil", "ciu", "ci_width", "cv")

  stopifnot(
    "jaccard_result must be a list with 'summary' and 'matrix_mean'" =
      is.list(jaccard_result) &&
      all(c("summary", "matrix_mean") %in% names(jaccard_result)),
    "stat must be one of the valid statistics" = stat %in% valid_stats
  )

  plot_df <- jaccard_result$summary %>%
    dplyr::mutate(stat_value = .data[[stat]])

  use_fixed_diagonal <- diagonal != "value"

  plot_df <- plot_df %>%
    dplyr::mutate(
      is_diagonal = unit_a == unit_b,
      fill_value  = dplyr::case_when(
        is_diagonal & use_fixed_diagonal ~ NA_real_,
        TRUE                             ~ stat_value
      ),
      label_value = dplyr::case_when(
        is_diagonal & use_fixed_diagonal ~ NA_real_,
        TRUE                             ~ stat_value
      )
    )

  unit_order <- rownames(jaccard_result$matrix_mean)
  plot_df <- plot_df %>%
    dplyr::mutate(
      unit_a = factor(unit_a, levels = unit_order),
      unit_b = factor(unit_b, levels = rev(unit_order))
    )

  plot_df_diag  <- plot_df %>% dplyr::filter(is_diagonal)
  plot_df_other <- plot_df %>% dplyr::filter(!is_diagonal)

  scale_limits <- if (stat %in% c("sd", "ci_width", "cv")) {
    c(0, max(plot_df$stat_value, na.rm = TRUE))
  } else {
    c(0, 1)
  }

  p <- ggplot2::ggplot(plot_df_other,
                       ggplot2::aes(x = unit_a, y = unit_b,
                                    fill = fill_value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5)

  if (use_fixed_diagonal) {
    p <- p +
      ggplot2::geom_tile(data        = plot_df_diag,
                         ggplot2::aes(x = unit_a, y = unit_b),
                         fill        = diagonal,
                         color       = "white",
                         linewidth   = 0.5,
                         inherit.aes = FALSE)
  }

  p <- p +
    ggplot2::scale_fill_gradient(
      limits   = scale_limits,
      low      = low_color,
      high     = high_color,
      na.value = "white",
      name     = paste0("Jaccard\n", stat)
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(size = 13, angle = 45,
                                              vjust = 1.1, hjust = 1,
                                              color = text_color),
      axis.text.y     = ggplot2::element_text(size = 13, color = text_color),
      legend.position = "right",
      legend.title    = ggplot2::element_text(size = 14),
      legend.text     = ggplot2::element_text(size = 12),
      panel.grid      = ggplot2::element_blank()
    )

  if (show_values) {
    plot_df_labels <- plot_df %>% dplyr::filter(!is.na(label_value))
    p <- p +
      ggplot2::geom_text(
        data        = plot_df_labels,
        ggplot2::aes(x = unit_a, y = unit_b,
                     label = round(label_value, digits)),
        color       = "black",
        size        = 4,
        inherit.aes = FALSE
      )
  }

  p
}


#' Summarise intra- and inter-group Jaccard from site-level pairs
#'
#' Takes site-level Jaccard output from \code{compute_jaccard()} and a named
#' vector mapping sites to groups, then classifies each site pair as within-
#' or between-group. Self-comparisons and duplicate pairs are removed.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta") and Thiago Bicudo
#'
#' @param jacc_sites Output from \code{compute_jaccard()} at site level
#'   (i.e., with \code{group_var = NULL}).
#' @param group_var Named character vector where names are site names and
#'   values are group labels (e.g., region).
#' @param dissimilarity Logical. If \code{TRUE}, converts Jaccard similarity
#'   to dissimilarity via \code{1 - J}. Credible interval bounds are inverted
#'   accordingly. Defaults to \code{FALSE}.
#'
#' @return A \code{tibble} with one row per unique site pair and columns:
#'   \code{pair}, \code{unit_a}, \code{unit_b}, \code{mean}, \code{median},
#'   \code{sd}, \code{cil}, \code{ciu}, \code{ci_width}, \code{cv},
#'   \code{group_a}, \code{group_b}, \code{comparison_type}
#'   (\code{"within"} or \code{"between"}).
#'
#' @seealso \code{\link{compute_jaccard}}, \code{\link{plot_jaccard_groups}}
#'
#' @export
summarise_jaccard_groups <- function(jacc_sites,
                                     group_var,
                                     dissimilarity = FALSE) {

  stopifnot(
    "jacc_sites must be a list with 'summary'" = is.list(jacc_sites) &&
      "summary" %in%
      names(jacc_sites),
    "group_var must be a named character vector" = is.character(group_var) &&
      !is.null(names(group_var)),
    "dissimilarity must be logical"              = is.logical(dissimilarity)
  )

  group_df <- tibble::tibble(
    site  = names(group_var),
    group = as.character(group_var)
  )

  pair_summary <- jacc_sites$summary %>%
    dplyr::filter(unit_a != unit_b) %>%
    dplyr::filter(unit_a < unit_b) %>%
    dplyr::left_join(group_df, by = c("unit_a" = "site")) %>%
    dplyr::rename(group_a = group) %>%
    dplyr::left_join(group_df, by = c("unit_b" = "site")) %>%
    dplyr::rename(group_b = group) %>%
    dplyr::mutate(
      comparison_type = dplyr::if_else(group_a == group_b, "within", "between")
    )

  if (dissimilarity) {
    pair_summary <- pair_summary %>%
      dplyr::mutate(
        mean    = 1 - mean,
        median  = 1 - median,
        cil_new = 1 - ciu,
        ciu_new = 1 - cil,
        cil     = cil_new,
        ciu     = ciu_new
      ) %>%
      dplyr::select(-cil_new, -ciu_new)
  }

  pair_summary
}


#' Dot plot of within- vs between-group Jaccard similarity
#'
#' Generates a strip/dot plot where each point is a site pair and large
#' colored points show the posterior mean with 95% credible interval for
#' within-group and between-group comparisons. Groups can be shown all at
#' once or filtered to a single focal group.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta") and Thiago Bicudo
#'
#' @param pair_summary Output from \code{summarise_jaccard_groups()}.
#' @param focal_group Character or \code{NULL}. If provided, shows only pairs
#'   involving this group on the y-axis (within + all between-group pairs).
#'   If \code{NULL} (default), all within-group and between-group combinations
#'   are shown.
#' @param dissimilarity Logical. If \code{TRUE}, adjusts the x-axis label to
#'   reflect the dissimilarity scale. Defaults to \code{FALSE}.
#' @param point_color_within Character. Hex color for within-group mean point
#'   and CI. Defaults to \code{"#EC0099"}.
#' @param point_color_between Character. Hex color for between-group mean point
#'   and CI. Defaults to \code{"#008AA1"}.
#' @param point_color_pairs Character. Hex color for individual site pair
#'   points. Defaults to \code{"#D99E73"}.
#' @param title Character or \code{NULL}. Plot title. Defaults to \code{NULL}.
#' @param text_color Character. Hex color for axis text. Defaults to
#'   \code{"#1B3644"}.
#'
#' @return A \code{ggplot2} object.
#'
#' @seealso \code{\link{summarise_jaccard_groups}}, \code{\link{compute_jaccard}}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar scale_color_manual labs theme_bw theme element_text element_line element_blank position_jitter
#'
#' @export
plot_jaccard_groups <- function(pair_summary,
                                focal_group         = NULL,
                                dissimilarity       = FALSE,
                                point_color_within  = "#EC0099",
                                point_color_between = "#008AA1",
                                point_color_pairs   = "#D99E73",
                                title               = NULL,
                                text_color          = "#1B3644") {

  stopifnot(
    "pair_summary must be a data frame" = is.data.frame(pair_summary),
    "required columns missing from pair_summary" = all(
      c("mean", "cil", "ciu", "comparison_type", "group_a", "group_b") %in%
        names(pair_summary)
    )
  )

  x_label <- if (dissimilarity) {
    "Jaccard dissimilarity (1 - J)"
  } else {
    "Jaccard similarity"
  }

  if (!is.null(focal_group)) {
    plot_df <- pair_summary %>%
      dplyr::filter(group_a == focal_group | group_b == focal_group) %>%
      dplyr::mutate(
        y_label = dplyr::if_else(
          comparison_type == "within",
          paste0(focal_group, " (within)"),
          dplyr::if_else(group_a == focal_group, group_b, group_a)
        )
      )
  } else {
    plot_df <- pair_summary %>%
      dplyr::mutate(
        y_label = dplyr::if_else(
          comparison_type == "within",
          paste0(group_a, " (within)"),
          paste0(group_a, " \u00d7 ", group_b)
        )
      )
  }

  group_means <- plot_df %>%
    dplyr::group_by(y_label, comparison_type) %>%
    dplyr::summarise(
      grand_mean = mean(mean, na.rm = TRUE),
      grand_cil  = mean(cil,  na.rm = TRUE),
      grand_ciu  = mean(ciu,  na.rm = TRUE),
      .groups    = "drop"
    )

  ggplot2::ggplot(plot_df, ggplot2::aes(x = mean, y = y_label)) +
    ggplot2::geom_point(
      color    = point_color_pairs,
      size     = 1.8,
      alpha    = 0.4,
      position = ggplot2::position_jitter(height = 0.1, width = 0)
    ) +
    ggplot2::geom_point(
      data        = group_means,
      ggplot2::aes(x = grand_mean, y = y_label, color = comparison_type),
      size        = 4.5,
      shape       = 16,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_errorbar(
      data        = group_means,
      ggplot2::aes(xmin  = grand_cil,
                   xmax  = grand_ciu,
                   y     = y_label,
                   color = comparison_type),
      width       = 0.2,
      linewidth   = 0.8,
      orientation = "y",
      inherit.aes = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c("within"  = point_color_within,
                 "between" = point_color_between),
      labels = c("within"  = "Within region",
                 "between" = "Between regions"),
      name   = NULL
    ) +
    ggplot2::labs(title = title, x = x_label, y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title         = ggplot2::element_text(size = 16, color = text_color),
      axis.text          = ggplot2::element_text(size = 13, color = text_color),
      panel.grid.major.y = ggplot2::element_line(color = "grey92"),
      panel.grid.major.x = ggplot2::element_line(color = "grey92"),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.position    = "bottom",
      legend.text        = ggplot2::element_text(size = 12)
    )
}
