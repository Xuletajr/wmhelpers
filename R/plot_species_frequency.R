#' Plot species frequency as a horizontal bar chart
#'
#' Creates a standardized horizontal bar chart to visualize species-level
#' frequencies such as number of detections or number of sites with detections.
#'
#' This helper is designed for ecological workflows where species metadata
#' (e.g. conservation status, endemism) are joined upstream.
#'
#' @param dat A data frame containing species-level summaries.
#' @param value Unquoted numeric variable indicating the frequency to plot.
#' @param species Unquoted variable indicating species names.
#' @param status Optional unquoted variable indicating conservation status
#'   (e.g. IUCN category).
#' @param flag Optional unquoted variable indicating a binary conservation flag
#'   (e.g. endemic or species of greatest conservation need).
#' @param status_ref Character string indicating the reference status
#'   (default = "LC").
#' @param flag_value Character string indicating the value used to activate
#'   the flag (default = "Yes").
#' @param flag_symbol Character symbol appended to labels when the flag is active.
#' @param xlab Character string for the x-axis label.
#' @param subtitle Character string for the plot subtitle.
#' @param breaks_x Optional numeric. Number of x-axis breaks.
#' @param limits_x Optional numeric vector of length 2 specifying x-axis limits.
#' @param config Optional list of project-level configuration parameters
#'   (commonly read from a YAML file in WildMon projects).
#'
#' @details
#' Label construction is controlled via \code{status_ref}, \code{flag_value},
#' and \code{flag_symbol}, allowing flexible representation of conservation
#' rules across projects.
#'
#' If \code{config} is provided, color settings are read from it; otherwise,
#' default ggplot2 styling is applied.
#'
#' @return A \code{ggplot} object.
#'
#' @author Jose W. Ribeiro Jr. ("Xuleta")
#'
#' @export
#'
plot_bar_species <- function(
    dat,
    value,
    species,
    status = NULL,
    flag = NULL,
    status_ref = "LC",
    flag_value = "Yes",
    flag_symbol = "*",
    xlab = NULL,
    subtitle = NULL,
    breaks_x = NULL,
    limits_x = NULL,
    config = NULL
) {

  # Resolve colors safely
  axis_color <- if (!is.null(config)) config$color2 else "black"
  bar_color  <- if (!is.null(config)) config$color7 else "grey70"

  # Build labels only if status is provided
  if (!is.null(status)) {
    dat <- dat %>%
      dplyr::mutate(
        label = dplyr::case_when(
          {{ status }} == status_ref &
            !is.null(flag) & {{ flag }} == flag_value ~ flag_symbol,

          {{ status }} != status_ref &
            !is.null(flag) & {{ flag }} == flag_value ~
            paste0({{ status }}, flag_symbol),

          TRUE ~ as.character({{ status }})
        )
      )
  } else {
    dat <- dat %>% dplyr::mutate(label = "")
  }

  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = {{ value }},
      y = reorder({{ species }}, {{ value }})
    )
  ) +
    ggplot2::geom_col(fill = bar_color) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      hjust = -0.2,
      size = 4
    ) +
    ggplot2::labs(
      x = xlab,
      y = NULL,
      subtitle = subtitle
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 14, color = axis_color),
      axis.text  = ggplot2::element_text(size = 12, color = axis_color),
      plot.subtitle = ggplot2::element_text(
        size = 13,
        color = axis_color,
        hjust = 1,
        face = "bold"
      ),
      legend.position = "none"
    )

  if (!is.null(breaks_x)) {
    p <- p +
      ggplot2::scale_x_continuous(
        n.breaks = breaks_x,
        limits = limits_x,
        labels = scales::label_number(accuracy = 1)
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

