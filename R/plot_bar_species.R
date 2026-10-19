#' Plot species frequency as a horizontal bar chart
#'
#' Creates a standardized horizontal bar chart to visualize species-level
#' frequencies (e.g., number of detections or number of sites).
#'
#' The function supports filtering to the top or bottom N species,
#' automatic formatting of scientific names (italic), and conservation
#' status labeling following WildMon standards. Full IUCN category names
#' are automatically converted to standard codes (e.g., "Least Concern"
#' to "LC").
#'
#' @param dat A data frame containing species-level summaries.
#' @param value Character string. Name of the numeric column used as the
#'   frequency metric (e.g., `"N_detections"` or `"N_sites"`).
#' @param species Character string. Name of the column containing species
#'   names (e.g., `"Species"` or `"Common_name"`).
#' @param status Optional character string. Name of the column containing
#'   conservation status (e.g., `"iucn_code"`). If provided, species with
#'   status matching `status_ref` are not labeled. Full IUCN category names
#'   are automatically converted to standard codes.
#' @param flag Optional character string. Name of a binary column used to
#'   highlight species of special interest (e.g., `"SGCN"`).
#' @param status_ref Character vector. Status category or categories that
#'   will not be displayed as labels (default = `c("LC", "NE")`).
#' @param flag_value Character string. Value in `flag` column used to
#'   activate the flag symbol (default = `"Yes"`).
#' @param flag_symbol Character string. Symbol appended to conservation
#'   labels when the flag is active (default = `"*"`).
#' @param top_n Optional numeric. If provided, selects only the top N
#'   (or bottom N if `bottom = TRUE`) species based on `value`.
#' @param bottom Logical. If `TRUE`, selects the least detected species
#'   instead of the most detected (default = `FALSE`).
#' @param bar_color Character string defining bar fill color
#'   (default = `"#00D973"`).
#' @param xlab Character string for the x-axis label.
#' @param subtitle Character string for the plot subtitle.
#' @param breaks_x Optional numeric. Number of breaks for the x-axis.
#' @param axis_text_size Numeric. Font size for axis labels and species names
#'   (default = `11`).
#' @param axis_title_size Numeric. Font size for axis titles
#'   (default = `12`).
#'
#' @details
#' If `species = "Species"`, names are formatted in italic.
#' If `species = "Common_name"`, names are displayed in plain text.
#'
#' Conservation labeling follows this logic:
#' - Full IUCN category names are automatically recoded to standard
#'   abbreviations (CR, EN, VU, NT, LC, DD, NE).
#' - Species with status matching any value in `status_ref` are not labeled.
#' - Species with a different status are labeled with their category code.
#' - If `flag` is provided and equals `flag_value`, the `flag_symbol`
#'   is appended to the label.
#'
#' @return A `ggplot2` object.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @export
plot_bar_species <- function(
    dat,
    value,
    species,
    status = NULL,
    flag = NULL,
    status_ref = c("LC", "NE"),
    flag_value = "Yes",
    flag_symbol = "*",
    top_n = NULL,
    bottom = FALSE,
    bar_color = "#00D973",
    xlab = NULL,
    subtitle = NULL,
    breaks_x = NULL,
    axis_text_size = 11,
    axis_title_size = 12
) {

  # ── Select top or bottom N species ───────────────────────────────────────────
  if (!is.null(top_n)) {
    if (bottom) {
      dat <- dat |>
        dplyr::arrange(.data[[value]]) |>
        dplyr::slice_head(n = top_n)
    } else {
      dat <- dat |>
        dplyr::arrange(dplyr::desc(.data[[value]])) |>
        dplyr::slice_head(n = top_n)
    }
  }

  # ── Recode full IUCN names to standard codes ─────────────────────────────────
  iucn_recode <- c(
    "Critically Endangered" = "CR",
    "Endangered"            = "EN",
    "Vulnerable"            = "VU",
    "Near Threatened"       = "NT",
    "Least Concern"         = "LC",
    "Data Deficient"        = "DD",
    "Not Evaluated"         = "NE"
  )

  if (!is.null(status) && status %in% names(dat)) {
    if (any(dat[[status]] %in% names(iucn_recode), na.rm = TRUE)) {
      dat[[status]] <- dplyr::recode(dat[[status]], !!!iucn_recode)
    }
  }

  # ── Build conservation label ──────────────────────────────────────────────────
  if (!is.null(status)) {

    if (!is.null(flag)) {

      dat$label <- dplyr::case_when(
        dat[[status]] %in% status_ref & dat[[flag]] == flag_value ~ flag_symbol,
        !dat[[status]] %in% status_ref & dat[[flag]] == flag_value ~
          paste0(dat[[status]], flag_symbol),
        !dat[[status]] %in% status_ref ~ dat[[status]],
        TRUE ~ ""
      )

    } else {

      dat$label <- ifelse(
        dat[[status]] %in% status_ref,
        "",
        dat[[status]]
      )
    }

  } else {
    dat$label <- ""
  }

  # ── Determine italic formatting for scientific names ──────────────────────────
  use_italic <- species == "Species"

  # ── Build plot ────────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = .data[[value]],
      y = reorder(.data[[species]], .data[[value]])
    )
  ) +
    ggplot2::geom_col(fill = bar_color) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      hjust = -0.2,
      size  = 4
    ) +
    ggplot2::labs(
      x        = xlab,
      y        = NULL,
      subtitle = subtitle
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(size = axis_text_size),
      axis.text.y     = ggplot2::element_text(
        size = axis_text_size,
        face = ifelse(use_italic, "italic", "plain")
      ),
      axis.title.x    = ggplot2::element_text(size = axis_title_size),
      legend.position = "none"
    )

  if (!is.null(breaks_x)) {
    p <- p +
      ggplot2::scale_x_continuous(
        n.breaks = breaks_x,
        labels   = scales::label_number(accuracy = 1)
      )
  } else {
    p <- p +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(c(0.01, 0.1)),
        labels = scales::label_number(accuracy = 1)
      )
  }

  return(p)
}


#' Plot detection summary for species
#'
#' High-level function that creates detection and sites plots for species-level
#' summaries. Produces paired bar charts showing number of detections and number
#' of sites for the top species.
#'
#' @param det_summary A data frame of detection summaries with at minimum
#'   columns `common_name`, status column (default `iucn_code`), `N_detections`,
#'   and `N_sites`.
#' @param top_n Integer. Number of top species to display (default = `10`).
#' @param color_det Character. Bar fill color for detections plot
#'   (default = `"#D99E73"`).
#' @param color_sites Character. Bar fill color for sites plot
#'   (default = `"#008AA1"`).
#' @param subtitle Character. Subtitle for the plots.
#' @param status_col Character. Name of the conservation status column
#'   (default = `"iucn_code"`).
#' @param axis_text_size Numeric. Font size for axis labels and species names
#'   (default = `11`).
#' @param axis_title_size Numeric. Font size for axis titles
#'   (default = `12`).
#'
#' @return A list with two `ggplot2` objects: `$detections` and `$sites`.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @export
plot_detection_summary <- function(det_summary,
                                   top_n            = 10,
                                   color_det        = "#D99E73",
                                   color_sites      = "#008AA1",
                                   subtitle         = NULL,
                                   status_col       = "iucn_code",
                                   axis_text_size   = 11,
                                   axis_title_size  = 12) {

  # ── Input validation ─────────────────────────────────────────────────────────
  stopifnot(
    "det_summary must be a data frame" = is.data.frame(det_summary),
    "top_n must be a positive integer" = is.numeric(top_n) && top_n > 0
  )

  required_cols <- c("common_name", status_col, "N_detections", "N_sites")
  missing_cols  <- setdiff(required_cols, names(det_summary))
  if (length(missing_cols) > 0) {
    stop("Missing columns in det_summary: ",
         paste(missing_cols, collapse = ", "))
  }

  # ── Detections plot ───────────────────────────────────────────────────────────
  p_det <- plot_bar_species(
    dat              = det_summary,
    value            = "N_detections",
    species          = "common_name",
    status           = status_col,
    top_n            = top_n,
    bar_color        = color_det,
    xlab             = "Number of detections",
    subtitle         = subtitle,
    axis_text_size   = axis_text_size,
    axis_title_size  = axis_title_size
  )

  # ── Sites plot ────────────────────────────────────────────────────────────────
  p_sites <- plot_bar_species(
    dat              = det_summary,
    value            = "N_sites",
    species          = "common_name",
    status           = status_col,
    top_n            = top_n,
    bar_color        = color_sites,
    xlab             = "Number of sites",
    subtitle         = subtitle,
    axis_text_size   = axis_text_size,
    axis_title_size  = axis_title_size
  )

  list(detections = p_det, sites = p_sites)
}


#' Export detection summary plots
#'
#' High-level pipeline function that iterates over a list of filter
#' specifications, applies each filter to the detection data, calls
#' \code{\link{summarise_detections}} and \code{\link{plot_detection_summary}},
#' and exports each plot as a separate PNG file.
#'
#' @param detections A data frame of detections containing at minimum
#'   `sites`, `species`, and any columns referenced in `plot_specs` filters
#'   (e.g., `method`, `season`, `community`, `landcover`).
#' @param spp_list A species metadata data frame passed to
#'   \code{\link{summarise_detections}}.
#' @param plot_specs A list of plot specifications. Each element must be a
#'   list with two elements:
#'   \describe{
#'     \item{\code{filters}}{A named list of column-value pairs used to filter
#'       `detections` before summarising. Use `NULL` for no filtering.}
#'     \item{\code{label}}{Character. A short descriptor used in the subtitle
#'       and output filename (e.g., `"audio_wet"`).}
#'   }
#' @param output_path Character. Directory where PNG files are saved.
#' @param project_name Character. Appended to output filenames.
#' @param dpi Numeric. Output resolution (default = `600`).
#' @param color_det Character. Bar fill color for detections plot
#'   (default = `"#D99E73"`).
#' @param color_sites Character. Bar fill color for sites plot
#'   (default = `"#008AA1"`).
#' @param top_n Integer. Number of top species per plot (default = `10`).
#' @param width Numeric. Output image width in cm (default = `18`).
#' @param height Numeric. Output image height in cm (default = `16`).
#' @param species_col Character. Species column name in `detections`
#'   (default = `"species"`).
#' @param common_name_col Character. Common name column name in `spp_list`
#'   (default = `"common_name"`).
#' @param iucn_col Character. IUCN status column name in `spp_list`
#'   (default = `"iucn_code"`).
#' @param axis_text_size Numeric. Font size for axis labels and species names
#'   (default = `11`).
#' @param axis_title_size Numeric. Font size for axis titles
#'   (default = `12`).
#'
#' @return Invisibly returns `NULL`. Called for its side effect of
#'   writing PNG files to `output_path`.
#'
#' @details
#' Each spec produces two separate PNG files:
#' \itemize{
#'   \item `species_detection_det_<label>_<project_name>.png`
#'   \item `species_detection_sites_<label>_<project_name>.png`
#' }
#'
#' If a filter combination results in zero rows, the iteration is skipped
#' with a warning and processing continues with the next specification.
#'
#' @seealso \code{\link{summarise_detections}},
#'   \code{\link{plot_detection_summary}}
#'
#' @importFrom purrr walk
#' @importFrom grDevices png dev.off
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @export
export_detection_plots <- function(detections,
                                   spp_list,
                                   plot_specs,
                                   output_path,
                                   project_name,
                                   dpi             = 600,
                                   color_det       = "#D99E73",
                                   color_sites     = "#008AA1",
                                   top_n           = 10,
                                   width           = 18,
                                   height          = 16,
                                   species_col     = "species",
                                   common_name_col = "common_name",
                                   iucn_col        = "iucn_code",
                                   axis_text_size  = 11,
                                   axis_title_size = 12) {

  # ── Input validation ─────────────────────────────────────────────────────────
  stopifnot(
    "detections must be a data frame"      = is.data.frame(detections),
    "spp_list must be a data frame"        = is.data.frame(spp_list),
    "plot_specs must be a list"            = is.list(plot_specs),
    "output_path must be a single string"  = is.character(output_path) &&
      length(output_path) == 1,
    "project_name must be a single string" = is.character(project_name) &&
      length(project_name) == 1,
    "dpi must be a positive number"        = is.numeric(dpi) && dpi > 0,
    "axis_text_size must be positive"      = is.numeric(axis_text_size) &&
      axis_text_size > 0,
    "axis_title_size must be positive"     = is.numeric(axis_title_size) &&
      axis_title_size > 0
  )

  if (!endsWith(output_path, "/")) output_path <- paste0(output_path, "/")

  # ── Validate plot_specs structure ─────────────────────────────────────────────
  for (i in seq_along(plot_specs)) {
    spec <- plot_specs[[i]]
    if (!all(c("filters", "label") %in% names(spec))) {
      stop("plot_specs[[", i, "]] must have 'filters' and 'label' elements.")
    }
  }

  # ── Iterate and export ────────────────────────────────────────────────────────
  purrr::walk(plot_specs, function(spec) {

    # ── Apply all filters dynamically ──────────────────────────────────────────
    dat <- detections

    if (!is.null(spec$filters) && length(spec$filters) > 0) {
      for (col in names(spec$filters)) {
        val <- spec$filters[[col]]
        if (!col %in% names(dat)) {
          warning("Filter column '", col, "' not found in detections. ",
                  "Skipping spec '", spec$label, "'.")
          return(invisible(NULL))
        }
        dat <- dplyr::filter(dat, .data[[col]] == val)
      }
    }

    if (nrow(dat) == 0) {
      warning("No rows remaining after filtering for spec '",
              spec$label, "'. Skipping.")
      return(invisible(NULL))
    }

    # ── Summarise ──────────────────────────────────────────────────────────────
    det_summary <- summarise_detections(
      detections_df   = dat,
      spp_list        = spp_list,
      species_col     = species_col,
      common_name_col = common_name_col,
      iucn_col        = iucn_col
    )

    # ── Build plots ────────────────────────────────────────────────────────────
    plots <- plot_detection_summary(
      det_summary     = det_summary,
      top_n           = top_n,
      color_det       = color_det,
      color_sites     = color_sites,
      subtitle        = spec$label,
      status_col      = "iucn_code",
      axis_text_size  = axis_text_size,
      axis_title_size = axis_title_size
    )

    # ── Export detections plot ─────────────────────────────────────────────────
    out_det <- paste0(
      output_path,
      "species_detection_det_", spec$label, "_", project_name, ".png"
    )
    grDevices::png(filename = out_det, width = width, height = height,
                   units = "cm", bg = "white", res = dpi)
    print(plots$detections)
    grDevices::dev.off()
    message("Exported: ", basename(out_det))

    # ── Export sites plot ──────────────────────────────────────────────────────
    out_sites <- paste0(
      output_path,
      "species_detection_sites_", spec$label, "_", project_name, ".png"
    )
    grDevices::png(filename = out_sites, width = width, height = height,
                   units = "cm", bg = "white", res = dpi)
    print(plots$sites)
    grDevices::dev.off()
    message("Exported: ", basename(out_sites))
  })

  invisible(NULL)
}
