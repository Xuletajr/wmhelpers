#' Export Occupancy Prediction Plots to PNG
#'
#' Saves individual PNG files for each occupancy prediction plot and
#' optionally creates a combined multiplot for a subset of variables.
#'
#' @param plot_list Named list of \code{ggplot} objects returned by
#'   \code{\link{plot_occ_predictions}}.
#' @param output_dir Character. Directory where plots will be saved.
#'   A trailing "/" will be added if missing.
#' @param project_name Character. String appended to all output filenames.
#' @param tag Optional character string inserted into all filenames before
#'   \code{.png} (e.g. \code{"audio"}, \code{"birds"}, \code{"amphibians"}).
#'   Defaults to \code{NULL} (no tag added).
#' @param type_levels Optional character vector (used only for filename tagging).
#' @param selected_vars Optional character vector of variables to include in
#'   the multiplot. If \code{NULL}, all plots are used.
#' @param width_single,height_single Numeric. Dimensions (cm) for individual plots.
#' @param width_multi,height_multi Numeric. Dimensions (cm) for multiplot.
#' @param dpi Integer. Resolution in dots per inch.
#'
#' @return Invisibly returns the multiplot object (if created), otherwise \code{NULL}.
#'
#' @description
#' Exports all plots in \code{plot_list} as individual PNG files and optionally
#' creates a combined multiplot using \code{patchwork}.
#'
#' File naming convention:
#' \code{occ_<variable>[_<type_tag>][_<tag>]_<project_name>.png}
#'
#' @author Thiago Bicudo
#'
#' @importFrom ggplot2 ggsave
#' @importFrom patchwork wrap_plots
#'
#' @export
export_occ_plots <- function(
    plot_list,
    output_dir,
    project_name,
    tag = NULL,
    type_levels = NULL,
    selected_vars = NULL,
    width_single = 24,
    height_single = 18,
    width_multi = 30,
    height_multi = 20,
    dpi = 600
) {

  # ── Input validation ───────────────────────────────────────────────
  stopifnot(
    "plot_list must be a list"         = is.list(plot_list),
    "output_dir must be character"     = is.character(output_dir),
    "project_name must be character"   = is.character(project_name)
  )

  if (!is.null(tag)) {
    stopifnot("tag must be a single character string" =
                is.character(tag) && length(tag) == 1)
  }

  # ── Ensure output directory ends with "/" ──────────────────────────
  if (!grepl("/$", output_dir)) {
    output_dir <- paste0(output_dir, "/")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # ── Build filename suffixes ────────────────────────────────────────
  type_tag <- if (!is.null(type_levels)) {
    paste0("_", paste(type_levels, collapse = "_"))
  } else {
    ""
  }

  custom_tag <- if (!is.null(tag)) paste0("_", tag) else ""

  # ── Save individual plots ──────────────────────────────────────────
  for (nm in names(plot_list)) {

    file_name <- paste0(
      output_dir,
      "occ_", nm, type_tag, custom_tag, "_", project_name, ".png"
    )

    ggplot2::ggsave(
      filename = file_name,
      plot     = plot_list[[nm]],
      width    = width_single,
      height   = height_single,
      units    = "cm",
      dpi      = dpi
    )
  }

  # ── Build multiplot (optional) ─────────────────────────────────────
  multi_plot <- NULL

  selected_plots <- if (!is.null(selected_vars)) {
    plot_list[
      sapply(names(plot_list), function(nm) {
        any(grepl(paste(selected_vars, collapse = "|"), nm))
      })
    ]
  } else {
    plot_list
  }

  if (length(selected_plots) > 1) {

    multi_plot <- patchwork::wrap_plots(selected_plots)

    multi_name <- paste0(
      output_dir,
      "occ_MULTIPLOT", type_tag, custom_tag, "_", project_name, ".png"
    )

    ggplot2::ggsave(
      filename = multi_name,
      plot     = multi_plot,
      width    = width_multi,
      height   = height_multi,
      units    = "cm",
      dpi      = dpi
    )
  }

  return(invisible(multi_plot))
}
