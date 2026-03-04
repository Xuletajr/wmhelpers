#' Plot species frequency as a horizontal bar chart
#'
#' Creates a standardized horizontal bar chart to visualize species-level
#' frequencies (e.g., number of detections or number of sites).
#'
#' The function supports filtering to the top or bottom N species,
#' automatic formatting of scientific names (italic), and conservation
#' status labeling following WildMon standards.
#'
#' @param dat A data frame containing species-level summaries.
#' @param value Character string. Name of the numeric column used as the
#'   frequency metric (e.g., `"N_detections"` or `"N_sites"`).
#' @param species Character string. Name of the column containing species
#'   names (e.g., `"Species"` or `"Common_name"`).
#' @param status Optional character string. Name of the column containing
#'   conservation status (e.g., `"IUCN"`). If provided, `"LC"` values
#'   are not displayed in the plot.
#' @param flag Optional character string. Name of a binary column used to
#'   highlight species of special interest (e.g., `"SGCN"`).
#' @param status_ref Character string. Reference status category that will
#'   not be displayed (default = `"LC"`).
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
#'
#' @details
#' If `species = "Species"`, names are formatted in italic.
#' If `species = "Common_name"`, names are displayed in plain text.
#'
#' Conservation labeling follows this logic:
#' - Species with status equal to `status_ref` (default `"LC"`)
#'   are not labeled unless flagged.
#' - Species with a different status are labeled with their category.
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
    status_ref = "LC",
    flag_value = "Yes",
    flag_symbol = "*",
    top_n = NULL,
    bottom = FALSE,
    bar_color = "#00D973",
    xlab = NULL,
    subtitle = NULL,
    breaks_x = NULL
) {

  # Select top or bottom N species
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

  # Build conservation label
  if (!is.null(status)) {

    if (!is.null(flag)) {

      dat$label <- dplyr::case_when(
        dat[[status]] == status_ref & dat[[flag]] == flag_value ~ flag_symbol,
        dat[[status]] != status_ref & dat[[flag]] == flag_value ~
          paste0(dat[[status]], flag_symbol),
        dat[[status]] != status_ref ~ dat[[status]],
        TRUE ~ ""
      )

    } else {

      dat$label <- ifelse(
        dat[[status]] == status_ref,
        "",
        dat[[status]]
      )
    }

  } else {
    dat$label <- ""
  }

  # Determine italic formatting for scientific names
  use_italic <- species == "Species"

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
      size = 4
    ) +
    ggplot2::labs(
      x = xlab,
      y = NULL,
      subtitle = subtitle
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(
        face = ifelse(use_italic, "italic", "plain")
      ),
      legend.position = "none"
    )

  if (!is.null(breaks_x)) {
    p <- p +
      ggplot2::scale_x_continuous(
        n.breaks = breaks_x,
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
