#' Export Site-Level Occupancy Probability Tables to CSV
#'
#' Writes CSV files for each occupancy probability summary statistic
#' (mean, SD, CRI 0.025, CRI 0.975) based on the output from
#' \code{\link{extract_psi_stats}}.
#'
#' @param psi_stats Named list returned by \code{\link{extract_psi_stats}}.
#'   Expected elements include matrices or data frames for:
#'   \code{mean}, \code{sd}, \code{criL}, and \code{criU}.
#' @param output_dir Character. Output directory. A trailing "/" is added
#'   if missing.
#' @param project_name Character. String appended to all filenames.
#'
#' @return Invisibly returns \code{NULL}. Function is called for its side effect.
#'
#' @description
#' Converts occupancy probability matrices into site × species tables and
#' exports one CSV per summary statistic. Each table includes a leading
#' \code{sites} column followed by species columns.
#'
#' @author Thiago Bicudo
#'
#' @importFrom utils write.csv
#'
#' @export
export_psi_tables <- function(psi_stats, output_dir, project_name) {

  # ── Input validation ───────────────────────────────────────────────
  stopifnot(
    "psi_stats must be a list" = is.list(psi_stats),
    "output_dir must be character" = is.character(output_dir),
    "project_name must be character" = is.character(project_name)
  )

  required_names <- c("mean", "sd", "criL", "criU")

  if (!all(required_names %in% names(psi_stats))) {
    stop("psi_stats must contain: mean, sd, criL, criU")
  }

  # ── Ensure output directory ends with "/" ──────────────────────────
  if (!grepl("/$", output_dir)) {
    output_dir <- paste0(output_dir, "/")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # ── Helper: format and export table ────────────────────────────────
  export_one <- function(mat, stat_name) {

    # Convert to data.frame
    df <- as.data.frame(mat)

    # Transpose: sites × species
    df_t <- as.data.frame(t(df))

    # Add sites column
    df_t$sites <- rownames(df_t)

    # Reorder columns
    df_t <- df_t[, c("sites", setdiff(names(df_t), "sites"))]

    # Clean rownames
    rownames(df_t) <- NULL

    # File name
    file_name <- paste0(
      output_dir,
      "psi_", stat_name, "_", project_name, ".csv"
    )

    utils::write.csv(df_t, file = file_name, row.names = FALSE)
  }

  # ── Export all statistics ──────────────────────────────────────────
  export_one(psi_stats$mean, "mean")
  export_one(psi_stats$sd,   "sd")
  export_one(psi_stats$criL, "cri025")
  export_one(psi_stats$criU, "cri975")

  return(invisible(NULL))
}
