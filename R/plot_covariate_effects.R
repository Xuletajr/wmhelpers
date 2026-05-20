#' Plot Species-Level Marginal Effects of Covariates
#'
#' Generates and saves species-level marginal effect plots for a given
#' covariate from a fitted multi-species occupancy model. Each plot shows
#' posterior mean estimates with 95\% credible intervals per species, coloured
#' by effect direction. Supports optional faceting by any grouping variable
#' and optional splitting of the data into two subsets based on a specific
#' value of the facet variable.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta") and Thiago Bicudo
#'
#' @param data A data frame of species-level parameter estimates joined with
#'   species trait data. Must contain columns \code{param_spp_code},
#'   \code{common_name}, \code{Mean}, \code{CRI_low}, \code{CRI_high}, and
#'   \code{Effect}. When \code{facet_var} or \code{taxa} are used, must also
#'   contain the corresponding grouping columns (e.g., \code{guild},
#'   \code{habitat}, \code{group}).
#' @param covariate Character. The covariate prefix to filter by, matching
#'   the beginning of \code{param_spp_code} before the \code{"-"} separator
#'   (e.g., \code{"days"}, \code{"canopyHeight"}, \code{"I(days^2)"}). Must
#'   not be \code{NULL} or empty.
#' @param model_type Character. Label used in the output filename to indicate
#'   whether this is a detection or occupancy model (e.g., \code{"detection"},
#'   \code{"occupancy"}).
#' @param xlabel_dict A named character vector mapping covariate names to
#'   human-readable x-axis labels. Must be defined per project since
#'   covariates vary across analyses. If the covariate is not found in the
#'   dictionary, the raw covariate name is used as the axis label.
#'
#'   Example:
#'   \preformatted{
#'   xlabel_dict <- c(
#'     "days"       = "Julian Day - linear term",
#'     "I(days^2)"  = "Julian Day - quadratic term",
#'     "numD"       = "Number of camera trap nights",
#'     "canopy"     = "Canopy height (m)"
#'   )
#'   }
#' @param title Character or \code{NULL}. Currently reserved for future use;
#'   plot titles are generated automatically from \code{taxa} and subset
#'   information. Defaults to \code{NULL}.
#' @param facet_var Character or \code{NULL}. Column name in \code{data} to
#'   use for faceting (e.g., \code{"guild"}, \code{"habitat"}, \code{"group"}).
#'   Use \code{NULL} to suppress faceting. Defaults to \code{NULL}.
#' @param split_by Character or \code{NULL}. A specific value of
#'   \code{facet_var} used to split the data into two subsets, generating two
#'   separate plots: one for species matching \code{split_by} and one for all
#'   others. For example, \code{split_by = "Invertivore"} with
#'   \code{facet_var = "guild"} produces one plot for invertivores and one
#'   for all other guilds. Use \code{NULL} to suppress splitting.
#'   Defaults to \code{NULL}.
#' @param taxa Character or \code{NULL}. Taxonomic group to filter by (e.g.,
#'   \code{"Aves"}, \code{"Mammalia"}). Use \code{"ALL"} or \code{NULL} to
#'   include all species. Defaults to \code{NULL}.
#' @param dodge_labels Logical. If \code{TRUE}, y-axis labels are staggered
#'   using \code{guide_axis(n.dodge = 2)} to reduce overlap. Defaults to
#'   \code{FALSE}.
#' @param text.size Numeric. Base font size for axis text. Must be a positive
#'   number.
#' @param filename_prefix Character. Prefix for the output filename.
#'   Defaults to \code{"species"}.
#' @param tag Character or \code{NULL}. Optional string appended to the
#'   output filename to identify the detection method or data source
#'   (e.g., \code{"audio"}, \code{"camtrap"}). Useful in integrated models
#'   where the same covariate exists across methods. Defaults to \code{NULL}.
#' @param output_dir Character. Directory where PNG files are saved. Defaults
#'   to \code{getwd()}.
#' @param project_name Character. String appended to output filenames to
#'   identify the project. Only the first 10 characters are used. Defaults
#'   to \code{""} (no suffix).
#' @param dpi Numeric. Resolution of output PNG files in dots per inch.
#'   Defaults to \code{600}.
#' @param color_neg Character. Hex colour for species with a negative effect.
#'   Defaults to \code{"#D99E73"}.
#' @param color_none Character. Hex colour for species with no significant
#'   effect. Defaults to \code{"#1B3644"}.
#' @param color_pos Character. Hex colour for species with a positive effect.
#'   Defaults to \code{"#008AA1"}.
#' @param color_strip Character. Hex colour for facet strip backgrounds.
#'   Defaults to \code{"#008AA1"}.
#' @param color_text Character. Hex colour for facet strip text.
#'   Defaults to \code{"#0F261F"}.
#'
#' @return A list of \code{ggplot2} objects, one per subset generated. Each
#'   plot is also printed to the active graphics device and saved as a PNG
#'   file in \code{output_dir}. The list is returned invisibly.
#'
#' @details
#' Output filenames are constructed from \code{filename_prefix},
#' \code{model_type}, the cleaned covariate name, optional facet and group
#' suffixes, \code{project_name} (first 10 characters), and \code{tag}.
#' Special characters in the covariate name (parentheses, carets, spaces,
#' slashes, colons) are removed or replaced with underscores.
#'
#' When \code{split_by} is provided, the suffix of the split subset value
#' is used directly (e.g., \code{"Invertivore"}) and the complementary
#' subset is named \code{"Other_<facet_var>"} (e.g., \code{"Other_guild"}).
#'
#' High-resolution output (30 cm x 55 cm at \code{2 * dpi}) is used when
#' \code{taxa} is \code{NULL} or \code{"ALL"}, or when \code{split_by} is
#' provided. Standard resolution (12 x 18 inches at \code{dpi}) is used
#' otherwise.
#'
#' Typical usage in an integrated model workflow:
#' \preformatted{
#' xlabel_dict <- c(
#'   "days"       = "Julian Day - linear term",
#'   "I(days^2)"  = "Julian Day - quadratic term",
#'   "numD"       = "Number of camera trap nights"
#' )
#'
#' # Audio - all species by guild, split invertivores from other guilds
#' plot_covariate_effects(
#'   data         = alpha_cov_aud,
#'   covariate    = "days",
#'   model_type   = "detection",
#'   xlabel_dict  = xlabel_dict,
#'   facet_var    = "guild",
#'   split_by     = "Invertivore",
#'   taxa         = "ALL",
#'   text.size    = 9,
#'   tag          = "audio",
#'   output_dir   = config$output_figures,
#'   project_name = config$project_name,
#'   dpi          = config$dpi
#' )
#'
#' # Camera trap - mammals, no split
#' plot_covariate_effects(
#'   data         = alpha_cov_cam,
#'   covariate    = "numD",
#'   model_type   = "detection",
#'   xlabel_dict  = xlabel_dict,
#'   facet_var    = NULL,
#'   split_by     = NULL,
#'   taxa         = "Mammalia",
#'   text.size    = 9,
#'   tag          = "camtrap",
#'   output_dir   = config$output_figures,
#'   project_name = config$project_name,
#'   dpi          = config$dpi
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_vline scale_colour_manual scale_y_discrete labs theme_classic theme element_text element_rect element_blank facet_grid guide_axis
#' @importFrom dplyr filter mutate
#' @importFrom stringr str_detect str_to_sentence str_escape
#' @importFrom scales alpha
#' @importFrom grDevices png dev.off
#' @importFrom stats as.formula
#'
#' @export
plot_covariate_effects <- function(data,
                                   covariate,
                                   model_type,
                                   xlabel_dict,
                                   title           = NULL,
                                   facet_var       = NULL,
                                   split_by        = NULL,
                                   taxa            = NULL,
                                   dodge_labels    = FALSE,
                                   text.size,
                                   filename_prefix = "species",
                                   tag             = NULL,
                                   output_dir      = getwd(),
                                   project_name    = "",
                                   dpi             = 600,
                                   color_neg       = "#D99E73",
                                   color_none      = "#1B3644",
                                   color_pos       = "#008AA1",
                                   color_strip     = "#008AA1",
                                   color_text      = "#0F261F") {

  # ── Input validation ───────────────────────────────────────────────────────
  .validate_plot_inputs(data, covariate, text.size)

  stopifnot(
    "model_type must be a single string"           = is.character(model_type) &&
      length(model_type) == 1,
    "xlabel_dict must be a named character vector" = is.character(xlabel_dict) &&
      !is.null(names(xlabel_dict)),
    "output_dir must be a single string"           = is.character(output_dir) &&
      length(output_dir) == 1,
    "project_name must be a character"             = is.character(project_name),
    "dpi must be a positive number"                = is.numeric(dpi) && dpi > 0,
    "color_neg must be a character string"         = is.character(color_neg),
    "color_none must be a character string"        = is.character(color_none),
    "color_pos must be a character string"         = is.character(color_pos),
    "color_strip must be a character string"       = is.character(color_strip),
    "color_text must be a character string"        = is.character(color_text)
  )

  # Normalise output_dir — add trailing slash if missing
  if (!endsWith(output_dir, "/")) output_dir <- paste0(output_dir, "/")

  # ── Resolve x-axis label from project dictionary ──────────────────────────
  x_axis_label <- if (covariate %in% names(xlabel_dict)) {
    xlabel_dict[[covariate]]
  } else {
    warning(
      "Covariate '", covariate, "' not found in xlabel_dict. ",
      "Using covariate name as axis label. ",
      "Add it to xlabel_dict to suppress this warning."
    )
    covariate
  }

  # ── Store original taxa for filename; NULL out "ALL" internally ────────────
  taxa_original <- taxa
  if (!is.null(taxa) && taxa == "ALL") taxa <- NULL

  # ── Filter and prepare data ────────────────────────────────────────────────
  filtered_data <- .prepare_plot_data(data, covariate, taxa)

  # ── Generate plot subsets ──────────────────────────────────────────────────
  plot_subsets <- .get_plot_subsets(filtered_data, facet_var, split_by)

  # ── Build, save, and collect plots ────────────────────────────────────────
  plots <- vector("list", length(plot_subsets))

  for (i in seq_along(plot_subsets)) {

    subset <- plot_subsets[[i]]

    plots[[i]] <- .create_and_save_plot(
      data_subset        = subset$data,
      suffix             = subset$suffix,
      covariate          = covariate,
      x_axis_label       = x_axis_label,
      facet_var          = facet_var,
      taxa_original      = taxa_original,
      model_type         = model_type,
      dodge_labels       = dodge_labels,
      base_text_size     = text.size,
      text_size_override = subset$text_size_override,
      title              = title,
      filename_prefix    = filename_prefix,
      tag                = tag,
      output_dir         = output_dir,
      project_name       = project_name,
      dpi                = dpi,
      color_neg          = color_neg,
      color_none         = color_none,
      color_pos          = color_pos,
      color_strip        = color_strip,
      color_text         = color_text
    )
  }

  for (p in plots) print(p)

  invisible(plots)
}


# =============================================================================
# Internal helper functions (not exported)
# =============================================================================

.validate_plot_inputs <- function(data, covariate, text.size) {

  if (is.null(data) || !is.data.frame(data)) {
    stop("data must be a non-NULL data frame.")
  }

  if (is.null(covariate) || !is.character(covariate) ||
      length(covariate) != 1 || nchar(covariate) == 0) {
    stop("covariate must be a single non-empty character string.")
  }

  if (!is.numeric(text.size) || length(text.size) != 1 || text.size <= 0) {
    stop("text.size must be a single positive number.")
  }

  required_cols <- c("param_spp_code", "common_name", "Mean",
                     "CRI_low", "CRI_high", "Effect")
  missing_cols  <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(
      "The following required columns are missing from data:\n",
      paste(missing_cols, collapse = ", ")
    )
  }
}


.prepare_plot_data <- function(data, covariate, taxa) {

  if (!is.null(taxa)) {
    data <- data %>% dplyr::filter(group == taxa)
  }

  covariate_pattern <- paste0("^", stringr::str_escape(covariate), "-")

  data %>%
    dplyr::filter(stringr::str_detect(param_spp_code, covariate_pattern)) %>%
    dplyr::mutate(common_name = stringr::str_to_sentence(common_name))
}


.get_plot_subsets <- function(data, facet_var, split_by) {

  # ── No split requested — return full data as single subset ─────────────────
  if (is.null(split_by) || is.null(facet_var)) {
    return(list(list(data = data, suffix = "", text_size_override = NULL)))
  }

  # ── Validate facet_var exists in data ──────────────────────────────────────
  if (!facet_var %in% names(data)) {
    stop(
      "facet_var '", facet_var, "' not found in data. ",
      "Cannot split by '", split_by, "'."
    )
  }

  # ── Validate split_by value exists in facet_var column ────────────────────
  available_vals <- unique(data[[facet_var]])

  if (!split_by %in% available_vals) {
    stop(
      "split_by value '", split_by, "' not found in column '",
      facet_var, "'.\n",
      "Available values: ",
      paste(sort(available_vals[!is.na(available_vals)]), collapse = ", ")
    )
  }

  # ── Build two subsets ─────────────────────────────────────────────────────
  suffix_in  <- split_by
  suffix_out <- paste0("Other_", facet_var)

  list(
    list(data               = dplyr::filter(data, .data[[facet_var]] == split_by),
         suffix             = suffix_in,
         text_size_override = NULL),
    list(data               = dplyr::filter(data, .data[[facet_var]] != split_by),
         suffix             = suffix_out,
         text_size_override = 14)
  )
}


.create_base_plot <- function(data, x_axis_label, dodge_labels, text_size,
                              plot_title,
                              color_neg, color_none, color_pos) {
  data %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(
      ggplot2::aes(x = Mean,
                   y = reorder(common_name, Mean),
                   colour = Effect),
      size = 4, shape = 16
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(y      = common_name,
                   xmin   = CRI_low,
                   xmax   = CRI_high,
                   colour = Effect),
      linewidth = 1.5, width = 0.3
    ) +
    ggplot2::scale_colour_manual(
      values = c("Negative"  = color_neg,
                 "No effect" = color_none,
                 "Positive"  = color_pos)
    ) +
    ggplot2::geom_vline(xintercept = 0, alpha = 0.75) +
    ggplot2::labs(
      x     = paste0("Parameter estimate (", x_axis_label, ")"),
      y     = NULL,
      title = plot_title
    ) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_discrete(
      guide = if (dodge_labels) {
        ggplot2::guide_axis(n.dodge = 2)
      } else {
        ggplot2::guide_axis()
      }
    ) +
    ggplot2::theme(
      axis.title.x      = ggplot2::element_text(size = 18),
      axis.text.x       = ggplot2::element_text(size = text_size),
      axis.text.y       = ggplot2::element_text(size   = text_size,
                                                hjust  = 1,
                                                colour = "black"),
      legend.background = ggplot2::element_rect(colour    = "black",
                                                linewidth = 0.5),
      legend.text       = ggplot2::element_text(size = 14),
      legend.position   = "top",
      legend.title      = ggplot2::element_blank(),
      plot.title        = ggplot2::element_text(size  = 20,
                                                hjust = 0.5,
                                                face  = "bold")
    )
}


.generate_plot_title <- function(title, taxa_original, suffix) {

  plot_title <- NULL

  if (!is.null(taxa_original) && taxa_original != "ALL") {
    plot_title <- paste("Taxonomic Group:", taxa_original)
  }

  if (!is.null(suffix) && suffix != "") {
    suffix_clean <- gsub("_", " ", suffix)
    plot_title   <- if (!is.null(plot_title)) {
      paste(plot_title, "-", suffix_clean)
    } else {
      suffix_clean
    }
  }

  plot_title
}


.add_facets <- function(plot, facet_var, color_strip, color_text) {

  if (!is.null(facet_var)) {
    plot <- plot +
      ggplot2::facet_grid(
        stats::as.formula(paste(facet_var, "~ .")),
        scales = "free_y",
        space  = "free_y"
      ) +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(
          fill = scales::alpha(color_strip, 0.6)
        ),
        strip.text.y = ggplot2::element_text(
          size   = 9,
          angle  = 0,
          colour = scales::alpha(color_text, 0.85),
          face   = "bold"
        )
      )
  }

  plot
}


.clean_covariate_name <- function(covariate) {
  covariate %>%
    gsub("\\\\", "", .) %>%
    gsub("\\(",  "", .) %>%
    gsub("\\)",  "", .) %>%
    gsub("\\^",  "", .) %>%
    gsub(" ",   "_", .) %>%
    gsub("/",   "_", .) %>%
    gsub(":",   "_", .)
}


.construct_filename <- function(model_type, covariate, facet_var,
                                taxa_original, suffix, filename_prefix,
                                tag, output_dir, project_name) {

  covariate_clean <- .clean_covariate_name(covariate)

  group_suffix <- if (!is.null(taxa_original) && taxa_original == "ALL") {
    "ALL"
  } else if (!is.null(taxa_original)) {
    taxa_original
  } else {
    ""
  }

  parts <- c(
    paste0(output_dir, filename_prefix, "_", model_type, "_", covariate_clean),
    if (!is.null(facet_var))                  facet_var                    else NULL,
    if (group_suffix != "")                   group_suffix                 else NULL,
    if (!is.null(suffix) && suffix != "")     suffix                       else NULL,
    if (nchar(project_name) > 0)              substr(project_name, 1, 10)  else NULL,
    if (!is.null(tag))                        tag                          else NULL
  )

  normalizePath(paste(parts, collapse = "_"), winslash = "/", mustWork = FALSE)
}


.should_use_high_quality <- function(taxa_original, split_by) {
  is.null(taxa_original) ||
    taxa_original == "ALL" ||
    !is.null(split_by)
}


.save_plot_to_file <- function(plot, filename, high_quality, dpi) {

  if (high_quality) {
    grDevices::png(paste0(filename, ".png"),
                   width     = 30 / 2.54,
                   height    = 55 / 2.54,
                   units     = "in",
                   pointsize = 6,
                   bg        = "white",
                   res       = dpi * 2)
  } else {
    grDevices::png(paste0(filename, ".png"),
                   width     = 12,
                   height    = 18,
                   units     = "in",
                   pointsize = 6,
                   bg        = "white",
                   res       = dpi)
  }

  print(plot)
  grDevices::dev.off()
}


.create_and_save_plot <- function(data_subset, suffix, covariate,
                                  x_axis_label, facet_var, taxa_original,
                                  model_type, dodge_labels, base_text_size,
                                  text_size_override = NULL, title = NULL,
                                  filename_prefix, tag,
                                  output_dir, project_name, dpi,
                                  color_neg, color_none, color_pos,
                                  color_strip, color_text) {

  text_size  <- if (!is.null(text_size_override)) text_size_override else base_text_size
  plot_title <- .generate_plot_title(title, taxa_original, suffix)

  plot <- .create_base_plot(data_subset, x_axis_label, dodge_labels,
                            text_size, plot_title,
                            color_neg, color_none, color_pos)

  plot <- .add_facets(plot, facet_var, color_strip, color_text)

  filename     <- .construct_filename(model_type, covariate, facet_var,
                                      taxa_original, suffix, filename_prefix,
                                      tag, output_dir, project_name)

  high_quality <- .should_use_high_quality(taxa_original, suffix)

  .save_plot_to_file(plot, filename, high_quality, dpi)

  message("Plot saved: ", basename(paste0(filename, ".png")))

  plot
}
