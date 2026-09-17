#' Export Site-Level Occupancy Probability Tables to CSV
#'
#' Writes CSV files for each occupancy probability summary statistic
#' (mean, SD, CRI 0.025, CRI 0.975) based on the output from
#' \code{\link{extract_psi_stats}}.
#'
#' @param psi_stats Named list returned by \code{\link{extract_psi_stats}}.
#'   Expected elements include matrices or data frames for:
#'   \code{mean}, \code{sd}, \code{cri025}, and \code{cri975}.
#' @param output_dir Character. Output directory. A trailing "/" is added
#'   if missing.
#' @param project_name Character. String appended to all filenames.
#'
#' @return Invisibly returns \code{NULL}. Function is called for its side effect.
#'
#' @description
#' Converts occupancy probability matrices into site x species tables and
#' exports one CSV per summary statistic. Each table includes a leading
#' \code{sites} column followed by species columns.
#'
#' @author Thiago Bicudo
#'
#' @importFrom utils write.csv
#'
#' @export
export_psi_tables <- function(psi_stats, output_dir, project_name) {
  stopifnot(
    `psi_stats must be a list` = is.list(psi_stats),
    `output_dir must be character` = is.character(output_dir),
    `project_name must be character` = is.character(project_name)
  )
  required_names <- c("mean", "sd", "cri025", "cri975")
  if (!all(required_names %in% names(psi_stats))) {
    stop("psi_stats must contain: mean, sd, cri025, cri975")
  }
  if (!grepl("/$", output_dir)) output_dir <- paste0(output_dir, "/")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  export_one <- function(mat, stat_name) {
    n_dims <- length(dim(mat))

    if (n_dims == 2) {
      # species x sites -> one row per site (unchanged from the original
      # single-season behaviour)
      df <- as.data.frame(t(mat))
      df$sites <- rownames(df)
      df <- df[, c("sites", setdiff(names(df), "sites"))]
      rownames(df) <- NULL
    } else if (n_dims == 3) {
      # species x sites x years -> one row per site-year, columns = species.
      # as.data.frame.table() gives each array dimension its OWN column (via
      # its dimnames), instead of relying on an implicit "site.year" string
      # concatenation that would need re-parsing downstream -- and would
      # break if a site name ever contained a period.
      long <- as.data.frame.table(mat, responseName = "value", stringsAsFactors = FALSE)
      names(long)[1:3] <- c("species", "sites", "year")
      df <- tidyr::pivot_wider(long, names_from = species, values_from = value) |>
        dplyr::relocate(sites, year)
    } else {
      stop("psi_stats element must be 2D (species x sites) or 3D ",
           "(species x sites x years); got ", n_dims, " dimensions.")
    }

    file_name <- paste0(output_dir, "psi_", stat_name, "_", project_name, ".csv")
    utils::write.csv(df, file = file_name, row.names = FALSE)
  }

  export_one(psi_stats$mean, "mean")
  export_one(psi_stats$sd, "sd")
  export_one(psi_stats$cri025, "cri025")
  export_one(psi_stats$cri975, "cri975")
  invisible(NULL)
}
