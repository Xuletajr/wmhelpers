# =============================================================================
# wmhelper — MSOM Covariate Effects & Marginal Effect Plots
#
# Functions in this file (in pipeline order):
#
#  COVARIATE LABELS
#   1. get_xlabel()              - human-readable x-axis label for a covariate
#   2. clean_covariate_name()    - sanitise covariate name for use in filenames
#
#  COVARIATE SUMMARIES
#   3. get_cov_info()            - CRI range + effect count summary per covariate
#   4. export_cov_info()         - write get_cov_info() results to a .txt file
#
#  MARGINAL EFFECT PLOTS
#   5. plot_covariate()          - main function: filter → plot → facet → save
#      Internal helpers (not exported):
#        .validate_plot_inputs()
#        .prepare_plot_data()
#        .get_plot_subsets()
#        .create_base_plot()
#        .generate_plot_title()
#        .add_facets()
#        .construct_filename()
#        .should_use_high_quality()
#        .save_plot_to_file()
#
#  FILE ORGANISATION
#   6. move_output_figures()     - move PNGs into structured subfolders
#
# =============================================================================
# WHAT TO KEEP AS PLAIN SCRIPT (not wrapped in functions):
#   - setwd / rm / gc / library calls
#   - config loading and model/data loading
#   - traits left_join
#   - defining det_covariates and occ_covariates vectors
#   - for() loops that call plot_covariate()
#   - dir.create() calls for new subfolders
#   - interactive browseURL() / print() checks
# =============================================================================


# =============================================================================
# 1. get_xlabel
# =============================================================================

#' Get Human-Readable X-Axis Label for a Covariate
#'
#' Returns a descriptive label for a model covariate name, suitable for use
#' as a plot axis title. Falls back to the raw covariate string if no label
#' is defined.
#'
#' @param covariate Character. The covariate name as it appears in
#'   \code{param_spp_code} (e.g. \code{"ndvi_mean"}, \code{"julian_std-"}).
#' @param custom_labels Optional named character vector of additional or
#'   override labels. Names are covariate strings, values are display labels.
#'   Merged with the built-in defaults, with \code{custom_labels} taking
#'   precedence. Defaults to \code{NULL}.
#'
#' @return A single character string.
#'
#' @examples
#' get_xlabel("ndvi_mean")
#' # "Normalized Difference Vegetation Index (NDVI)"
#'
#' get_xlabel("my_var", custom_labels = c(my_var = "My Custom Variable"))
#' # "My Custom Variable"
#'
#' @export
get_xlabel <- function(covariate, custom_labels = NULL) {

  default_labels <- c(
    "julian_std-"          = "Julian Day \u2013 linear term",
    "I\\(julian_std\\^2\\)"= "Julian Day \u2013 quadratic term",
    "canopy_mea"           = "Canopy Height (mean)",
    "canopy_std"           = "Canopy Height (standard deviation)",
    "mndwi_mean"           = "Modified Normalized Difference Water Index (MNDWI)",
    "ndvi_mean"            = "Normalized Difference Vegetation Index (NDVI)",
    "footprint"            = "Human Footprint Index",
    "bio7"                 = "Annual Temperature Range (Bio7)",
    "bio8"                 = "Mean Temperature of Wettest Quarter (Bio8)",
    "precipitation"        = "Average Daily Precipitation",
    "dist_water"           = "Distance to Water (natural log of distance)"
  )

  # Merge: custom_labels override defaults
  all_labels <- if (!is.null(custom_labels)) {
    c(custom_labels, default_labels[!names(default_labels) %in% names(custom_labels)])
  } else {
    default_labels
  }

  label <- all_labels[covariate]

  if (is.na(label) || is.null(label)) covariate else unname(label)
}


# =============================================================================
# 2. clean_covariate_name
# =============================================================================

#' Sanitise a Covariate Name for Use in Filenames
#'
#' Removes or replaces characters that are illegal or problematic in file
#' paths (\code{(}, \code{)}, \code{^}, \code{\\}, spaces, \code{/}).
#'
#' @param covariate Character. Raw covariate name.
#'
#' @return A sanitised character string safe for use in file paths.
#'
#' @examples
#' clean_covariate_name("I(julian_std^2)")
#' # "Ijulian_std2"
#'
#' @export
clean_covariate_name <- function(covariate) {
  covariate %>%
    gsub("\\\\", "", .) %>%   # backslashes
    gsub("\\(",  "", .) %>%   # opening parenthesis
    gsub("\\)",  "", .) %>%   # closing parenthesis
    gsub("\\^",  "", .) %>%   # caret
    gsub(" ",   "_", .) %>%   # spaces → underscores
    gsub("/",   "_", .)       # forward slash → underscore
}


# =============================================================================
# 3. get_cov_info
# =============================================================================

#' Summarise Covariate Effect Range and Direction
#'
#' Filters a parameter table to a specific covariate, computes the overall
#' CRI range (min of \code{CRI_low}, max of \code{CRI_high}), and returns
#' a count of effect directions by a chosen grouping variable.
#'
#' @param data Data frame containing at minimum columns \code{param_spp_code},
#'   \code{CRI_low}, \code{CRI_high}, and \code{Effect}, plus whatever column
#'   is named in \code{group_var}.
#' @param covariate Character. Regex pattern used to filter
#'   \code{param_spp_code} (e.g. \code{"ndvi_mean"} or
#'   \code{"I\\\\(julian_std\\\\^2\\\\)"}).
#' @param group_var Character. Column name to group the effect-direction
#'   summary by. Defaults to \code{"guild"}.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{\code{min}}{Minimum value of \code{CRI_low} across matched rows.}
#'   \item{\code{max}}{Maximum value of \code{CRI_high} across matched rows.}
#'   \item{\code{summary}}{Grouped tibble with effect counts.}
#' }
#'
#' @importFrom dplyr filter group_by summarise n
#' @importFrom stringr str_detect
#'
#' @export
get_cov_info <- function(data, covariate, group_var = "guild") {

  stopifnot(
    "data must be a data frame"          = is.data.frame(data),
    "covariate must be a character string" = is.character(covariate),
    "group_var must be a column in data" = group_var %in% names(data)
  )

  covariate_data <- data %>%
    dplyr::filter(stringr::str_detect(param_spp_code, covariate)) %>%
    dplyr::mutate(common_name = stringr::str_to_sentence(common_name))

  if (nrow(covariate_data) == 0) {
    warning("No rows matched covariate pattern: ", covariate)
    return(list(min = NA_real_, max = NA_real_, summary = NULL))
  }

  summary_data <- covariate_data %>%
    dplyr::group_by(.data[[group_var]], Effect) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  list(
    min     = min(covariate_data$CRI_low,  na.rm = TRUE),
    max     = max(covariate_data$CRI_high, na.rm = TRUE),
    summary = summary_data
  )
}


# =============================================================================
# 4. export_cov_info
# =============================================================================

#' Export Covariate Effect Summaries to a Text File
#'
#' Iterates over a vector of covariates, calls \code{\link{get_cov_info}}
#' for each, and writes the results to a plain-text file via \code{sink()}.
#'
#' @param data Data frame passed to \code{\link{get_cov_info}}.
#' @param covariates Character vector of covariate patterns to summarise.
#' @param output_path Character. Full path to the output \code{.txt} file.
#' @param model_type Character. Label printed in the file header, e.g.
#'   \code{"detection"} or \code{"occupancy"}.
#' @param group_var Character. Grouping variable for the summary.
#'   Defaults to \code{"guild"}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect.
#'
#' @export
export_cov_info <- function(data,
                             covariates,
                             output_path,
                             model_type = "occupancy",
                             group_var  = "guild") {

  stopifnot(
    "data must be a data frame"        = is.data.frame(data),
    "covariates must be a character vector" = is.character(covariates),
    "output_path must be a string"     = is.character(output_path)
  )

  sink(file = output_path)

  cat("\n",
      "########################################################################\n",
      "########  Marginal Effects of environmental variables  #################\n",
      "########################################################################\n",
      "Model type:", toupper(model_type), "\n\n\n")

  for (covariate in covariates) {

    cat("Marginal Effects of:", covariate, "\n")

    info <- get_cov_info(data, covariate, group_var)

    cat("Min CRI_low :", round(info$min, 4), "\n")
    cat("Max CRI_high:", round(info$max, 4), "\n")
    cat("\nSummary by", group_var, ":\n")
    print(info$summary, n = Inf)
    cat("------------------------------------------------------------------------\n\n")
  }

  sink()
  message("Covariate summary saved: ", basename(output_path))
  invisible(NULL)
}


# =============================================================================
# 5. plot_covariate   (+ internal helpers prefixed with .)
# =============================================================================

# ── Internal helpers ──────────────────────────────────────────────────────────

.validate_plot_inputs <- function(data, covariate, text_size) {
  if (is.null(data) || nrow(data) == 0)
    stop("'data' cannot be NULL or empty.", call. = FALSE)
  if (is.null(covariate) || covariate == "")
    stop("'covariate' must be specified.", call. = FALSE)
  if (!is.numeric(text_size) || text_size <= 0)
    stop("'text_size' must be a positive number.", call. = FALSE)
}

.prepare_plot_data <- function(data, covariate, taxa) {

  if (!is.null(taxa)) {
    data <- data %>% dplyr::filter(group == taxa)
  }

  data %>%
    dplyr::filter(stringr::str_detect(param_spp_code, covariate)) %>%
    dplyr::mutate(common_name = stringr::str_to_sentence(common_name))
}

.get_plot_subsets <- function(data, taxa, facet_var) {

  # No taxa filter → single subset
  if (is.null(taxa)) {
    return(list(list(data = data, suffix = "", text_size_override = NULL)))
  }

  # Birds with guild or habitat faceting → split into two subsets to keep
  # plot heights manageable
  if (taxa == "Aves" && !is.null(facet_var)) {

    if (facet_var == "guild") {
      return(list(
        list(data              = data %>% dplyr::filter(guild == "Invertivore"),
             suffix            = "Invertivore",
             text_size_override = NULL),
        list(data              = data %>% dplyr::filter(guild != "Invertivore"),
             suffix            = "Other_Guilds",
             text_size_override = 14)
      ))
    }

    if (facet_var == "habitat") {
      return(list(
        list(data              = data %>% dplyr::filter(habitat == "Forest"),
             suffix            = "Forest_Species",
             text_size_override = NULL),
        list(data              = data %>% dplyr::filter(habitat != "Forest"),
             suffix            = "Other_Species",
             text_size_override = 14)
      ))
    }
  }

  list(list(data = data, suffix = "", text_size_override = NULL))
}

.generate_plot_title <- function(taxa_original, suffix) {

  parts <- c(
    if (!is.null(taxa_original) && taxa_original != "ALL")
      paste("Taxonomic Group:", taxa_original),
    if (!is.null(suffix) && suffix != "")
      gsub("_", " ", suffix)
  )

  if (length(parts) == 0) NULL else paste(parts, collapse = " \u2013 ")
}

.create_base_plot <- function(data,
                               x_axis_label,
                               dodge_labels,
                               text_size,
                               plot_title,
                               color_negative,
                               color_noeffect,
                               color_positive) {

  data %>%
    ggplot2::ggplot() +
    ggplot2::geom_vline(
      xintercept = 0,
      colour     = "grey40",
      linewidth  = 0.6,
      alpha      = 0.8
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(y    = common_name,
                   xmin = CRI_low,
                   xmax = CRI_high,
                   colour = Effect),
      linewidth = 1.2,
      width     = 0.3
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x      = Mean,
                   y      = reorder(common_name, Mean),
                   colour = Effect),
      size  = 3.5,
      shape = 16
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "Negative"  = color_negative,
        "No effect" = color_noeffect,
        "Positive"  = color_positive
      )
    ) +
    ggplot2::labs(
      x     = paste0("Parameter estimate (", x_axis_label, ")"),
      y     = NULL,
      title = plot_title
    ) +
    ggplot2::scale_y_discrete(
      guide = if (dodge_labels) ggplot2::guide_axis(n.dodge = 2)
              else               ggplot2::guide_axis()
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.title.x     = ggplot2::element_text(size = 20),
      axis.text.x      = ggplot2::element_text(size = text_size),
      axis.text.y      = ggplot2::element_text(size  = text_size,
                                               hjust = 1,
                                               colour = "black"),
      legend.background = ggplot2::element_rect(colour    = "black",
                                                linewidth = 0.5),
      legend.text       = ggplot2::element_text(size = 14),
      legend.position   = "top",
      legend.title      = ggplot2::element_blank(),
      plot.title        = ggplot2::element_text(size  = 16,
                                                hjust = 0.5,
                                                face  = "bold")
    )
}

.add_facets <- function(plot, facet_var, color_strip_fill, color_strip_text) {

  if (!is.null(facet_var)) {
    plot <- plot +
      ggplot2::facet_grid(
        stats::as.formula(paste(facet_var, "~ .")),
        scales = "free_y",
        space  = "free_y"
      ) +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(
          fill = scales::alpha(color_strip_fill, 0.6)
        ),
        strip.text.y     = ggplot2::element_text(
          size   = 9,
          angle  = 0,
          colour = scales::alpha(color_strip_text, 0.85),
          face   = "bold"
        )
      )
  }

  plot
}

.construct_filename <- function(output_dir,
                                 model_type,
                                 covariate,
                                 facet_var,
                                 taxa_original,
                                 suffix,
                                 project_name,
                                 filename_prefix = "species") {

  covariate_clean <- clean_covariate_name(covariate)

  # Include facet_var in filename only when it wasn't used to split subsets
  bird_split_suffixes <- c("Invertivore", "Other_Guilds",
                           "Forest_Species", "Other_Species")
  facet_suffix <- if (!is.null(facet_var) &&
                       !suffix %in% bird_split_suffixes) facet_var else ""

  group_suffix <- if (!is.null(taxa_original) && taxa_original == "ALL") "ALL"
                  else if (!is.null(taxa_original))                       taxa_original
                  else                                                     ""

  parts <- c(
    paste0(output_dir, filename_prefix, "_", model_type, "_", covariate_clean),
    if (facet_suffix != "") facet_suffix,
    if (group_suffix != "") group_suffix,
    if (suffix       != "") suffix,
    project_name
  )

  normalizePath(paste(parts, collapse = "_"),
                winslash  = "/",
                mustWork  = FALSE)
}

.should_use_high_quality <- function(taxa, facet_var) {
  is.null(taxa) ||
    (!is.null(facet_var) && facet_var %in% c("habitat", "guild"))
}

.save_plot_to_file <- function(plot, filename, high_quality = FALSE, dpi = 300) {

  dims <- if (high_quality) {
    list(width = 30 / 2.54, height = 50 / 2.54, res = dpi * 2)
  } else {
    list(width = 12, height = 18, res = dpi)
  }

  png(paste0(filename, ".png"),
      width     = dims$width,
      height    = dims$height,
      units     = "in",
      pointsize = 6,
      bg        = "white",
      res       = dims$res)
  print(plot)
  dev.off()
}


# ── Main exported function ────────────────────────────────────────────────────

#' Plot Marginal Effects of a Covariate on Detection or Occupancy
#'
#' Filters a species-level parameter table to a specific covariate, builds
#' a dot-and-errorbar plot coloured by effect direction, optionally facets
#' by a trait variable, and saves the output to PNG. Handles automatic
#' splitting of large bird datasets into manageable subplots.
#'
#' @param data Data frame of species-level parameter estimates, typically
#'   produced by \code{\link{build_param_table}} joined with a species traits
#'   table. Must contain at minimum columns \code{param_spp_code},
#'   \code{common_name}, \code{Mean}, \code{CRI_low}, \code{CRI_high},
#'   \code{Effect}, and whatever column is named in \code{facet_var}.
#' @param covariate Character. Regex pattern matching the covariate in
#'   \code{param_spp_code} (e.g. \code{"ndvi_mean"},
#'   \code{"I\\\\(julian_std\\\\^2\\\\)"}).
#' @param output_dir Character. Directory for output PNG files.
#' @param project_name Character. Appended to all output filenames.
#' @param model_type Character. Used in the output filename and title
#'   (\code{"detection"} or \code{"occupancy"}).
#' @param title Character. Currently unused (reserved for future use);
#'   titles are generated automatically from \code{taxa} and \code{suffix}.
#' @param facet_var Character or \code{NULL}. Column in \code{data} to use
#'   for \code{facet_grid} rows (e.g. \code{"habitat"}, \code{"guild"}).
#'   Pass \code{NULL} to suppress faceting.
#' @param taxa Character or \code{NULL}. Value of the \code{group} column
#'   used to subset species. Pass \code{"ALL"} or \code{NULL} to include
#'   all species.
#' @param text_size Numeric. Base font size for axis tick labels.
#' @param dodge_labels Logical. If \code{TRUE}, y-axis labels are staggered
#'   to reduce overlap. Defaults to \code{FALSE}.
#' @param filename_prefix Character. Prefix for the output filename.
#'   Defaults to \code{"species"}. Use \code{"significant_species"} to
#'   distinguish filtered plots.
#' @param color_negative,color_noeffect,color_positive Characters.
#'   Colours for negative, neutral, and positive effects.
#' @param color_strip_fill,color_strip_text Characters. Colours for facet
#'   strip background and text.
#' @param dpi Integer. Resolution for PNG output. Defaults to \code{300}.
#' @param custom_labels Optional named character vector of x-axis label
#'   overrides passed to \code{\link{get_xlabel}}.
#'
#' @return Invisibly returns a list of \code{ggplot} objects (one per subset).
#'
#' @importFrom dplyr filter mutate
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_vline
#'   scale_colour_manual labs scale_y_discrete facet_grid theme_classic theme
#'   element_text element_rect element_blank guide_axis
#' @importFrom stringr str_detect str_to_sentence
#' @importFrom scales alpha
#'
#' @export
plot_covariate <- function(data,
                            covariate,
                            output_dir,
                            project_name,
                            model_type      = "occupancy",
                            title           = NULL,
                            facet_var       = NULL,
                            taxa            = NULL,
                            text_size       = 18,
                            dodge_labels    = FALSE,
                            filename_prefix = "species",
                            color_negative  = "#D95F02",
                            color_noeffect  = "#7570B3",
                            color_positive  = "#1B9E77",
                            color_strip_fill = "#CCCCCC",
                            color_strip_text = "#1B3A5C",
                            dpi             = 300,
                            custom_labels   = NULL) {

  # ── Validation ───────────────────────────────────────────────────────────────
  .validate_plot_inputs(data, covariate, text_size)
  if (!endsWith(output_dir, "/")) output_dir <- paste0(output_dir, "/")

  # ── Store original taxa for filename/title; treat "ALL" as no filter ─────────
  taxa_original <- taxa
  if (!is.null(taxa) && taxa == "ALL") taxa <- NULL

  # ── Filter and prepare data ───────────────────────────────────────────────────
  filtered_data <- .prepare_plot_data(data, covariate, taxa)

  if (nrow(filtered_data) == 0) {
    warning("No data remaining after filtering for covariate '", covariate,
            "' and taxa '", taxa_original, "'. Skipping.", call. = FALSE)
    return(invisible(NULL))
  }

  # ── X-axis label ──────────────────────────────────────────────────────────────
  x_axis_label <- get_xlabel(covariate, custom_labels)

  # ── Split into subsets (handles large bird datasets) ──────────────────────────
  subsets <- .get_plot_subsets(filtered_data, taxa, facet_var)

  # ── Build, save, and collect plots ───────────────────────────────────────────
  plots <- lapply(subsets, function(s) {

    ts    <- if (!is.null(s$text_size_override)) s$text_size_override else text_size
    ptitle <- .generate_plot_title(taxa_original, s$suffix)

    p <- .create_base_plot(
      data          = s$data,
      x_axis_label  = x_axis_label,
      dodge_labels  = dodge_labels,
      text_size     = ts,
      plot_title    = ptitle,
      color_negative  = color_negative,
      color_noeffect  = color_noeffect,
      color_positive  = color_positive
    )

    p <- .add_facets(p, facet_var, color_strip_fill, color_strip_text)

    fname <- .construct_filename(
      output_dir      = output_dir,
      model_type      = model_type,
      covariate       = covariate,
      facet_var       = facet_var,
      taxa_original   = taxa_original,
      suffix          = s$suffix,
      project_name    = project_name,
      filename_prefix = filename_prefix
    )

    hq <- .should_use_high_quality(taxa, facet_var)
    .save_plot_to_file(p, fname, hq, dpi)

    cat("Saved:", basename(fname), ".png\n", sep = "")
    p
  })

  # Print to screen
  for (p in plots) print(p)

  invisible(plots)
}


# =============================================================================
# 6. move_output_figures
# =============================================================================

#' Move Model Output Figures into Structured Subfolders
#'
#' Scans a base output directory for PNG files matching naming conventions
#' produced by \code{\link{plot_covariate}} and moves them into
#' \code{1.species/} and \code{2.community/} subfolders. Optionally moves
#' taxon-specific files into further sub-subfolders.
#'
#' @param base_figures_dir Character. Base figures directory (e.g.
#'   \code{config$output_figures}).
#' @param taxa_subdirs Optional named list mapping regex patterns to target
#'   sub-subdirectory names within \code{1.species/} and \code{2.community/}.
#'   Example:
#'   \code{list("_Birds_" = "birds", "_Frogs_" = "frogs")}.
#'   If \code{NULL} (default), taxon-level subdirectory moves are skipped.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect.
#'
#' @details
#' The function first closes any open graphics devices to ensure all PNG
#' files have been fully written before moving. Directory structure created:
#' \preformatted{
#'   Occupancy_Models/
#'     1.species/
#'       birds/      (if taxa_subdirs provided)
#'       frogs/      (if taxa_subdirs provided)
#'     2.community/
#'       birds/      (if taxa_subdirs provided)
#'       frogs/      (if taxa_subdirs provided)
#' }
#'
#' @importFrom fs dir_ls file_move
#'
#' @export
move_output_figures <- function(base_figures_dir, taxa_subdirs = NULL) {

  stopifnot(
    "base_figures_dir must be a character string" =
      is.character(base_figures_dir)
  )

  if (!endsWith(base_figures_dir, "/"))
    base_figures_dir <- paste0(base_figures_dir, "/")

  occ_dir     <- paste0(base_figures_dir, "Occupancy_Models/")
  species_dir <- paste0(occ_dir, "1.species/")
  comm_dir    <- paste0(occ_dir, "2.community/")

  # Close any open devices before moving
  while (!is.null(dev.list())) dev.off()

  to_norm <- function(path) normalizePath(path, winslash = "/", mustWork = FALSE)

  # ── Move species-level detection and occupancy figures ────────────────────────
  fs::dir_ls(occ_dir, regexp = "species_detection_.*\\.png$") %>%
    fs::file_move(to_norm(species_dir))

  fs::dir_ls(occ_dir, regexp = "species_occupancy_.*\\.png$") %>%
    fs::file_move(to_norm(species_dir))

  fs::dir_ls(occ_dir, regexp = "significant_species_.*\\.png$") %>%
    fs::file_move(to_norm(species_dir))

  # ── Move community figures ────────────────────────────────────────────────────
  fs::dir_ls(occ_dir, regexp = "community_.*\\.png$") %>%
    fs::file_move(to_norm(comm_dir))

  message("Species and community figures moved.")

  # ── Optionally move taxon-specific figures into sub-subfolders ───────────────
  if (!is.null(taxa_subdirs)) {

    for (pattern in names(taxa_subdirs)) {

      subdir_name <- taxa_subdirs[[pattern]]

      # Create sub-subfolders if they do not exist
      sp_subdir   <- to_norm(paste0(species_dir, subdir_name, "/"))
      comm_subdir <- to_norm(paste0(comm_dir,    subdir_name, "/"))
      dir.create(sp_subdir,   recursive = TRUE, showWarnings = FALSE)
      dir.create(comm_subdir, recursive = TRUE, showWarnings = FALSE)

      # Move matching files
      sp_matches <- fs::dir_ls(species_dir, regexp = paste0(pattern, ".*\\.png$"))
      if (length(sp_matches) > 0) fs::file_move(sp_matches, sp_subdir)

      comm_matches <- fs::dir_ls(comm_dir, regexp = paste0(pattern, ".*\\.png$"))
      if (length(comm_matches) > 0) fs::file_move(comm_matches, comm_subdir)

      message("Moved '", pattern, "' figures to: ", subdir_name)
    }
  }

  cat("All figures organised successfully.\n")
  invisible(NULL)
}
