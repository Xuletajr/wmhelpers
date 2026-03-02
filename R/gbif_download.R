# =============================================================================
# === FUNCTIONS. NO NEED TO EDIT IT ===========================================
# =============================================================================

# 2. Download Occurrence Data from GBIF -------------------------------------

# INTERNAL HELPER - Read species file (CSV or Excel) -------------------------
#' Detects file extension and reads accordingly. Used internally by
#' process_species_input and debug_name_matching.
#'
#' @param file_path Full path to .csv, .xls, or .xlsx file
#' @param sheet Sheet name (only used for Excel files, ignored for CSV)
#' @return A tibble with the file contents
read_species_file <- function(file_path, sheet = NULL) {
  
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  file_ext <- tolower(tools::file_ext(file_path))
  message("Detected file type: .", file_ext)
  
  species_data <- tryCatch({
    if (file_ext == "csv") {
      message("Reading CSV file...")
      readr::read_csv(file_path, show_col_types = FALSE)
      
    } else if (file_ext %in% c("xlsx", "xls")) {
      message("Reading Excel file...")
      readxl::read_xlsx(file_path, sheet = sheet)
      
    } else {
      stop("Unsupported file type: '.", file_ext,
           "'. Please provide a .csv, .xlsx, or .xls file.")
    }
  }, error = function(e) {
    stop("Error reading file: ", e$message)
  })
  
  return(species_data)
}

####
# 2.1. MAIN FUNCTION - Download GBIF Occurrence Data ------------------------
#' This function downloads species occurrence data from GBIF database.
#' It can handle both individual species lists and broad taxonomic groups.
#' All downloads include only records with geographic coordinates.
#'
#' @param species_input Character vector of species names OR "spreadsheet" to 
#'        read from a file (.csv, .xlsx, or .xls)
#' @param taxon_group Single taxonomic group name (e.g., "Aves", "Mammalia", "Amphibia"), 
#'        it prioritizes species_input when both parameters are provided
#' @param rank Taxonomic rank for taxon_group (e.g., "class", "order", "family"). 
#'        If NULL, function will try without rank first
#' @param country_code ISO 2-letter country code to filter results (e.g., "PR", "VI", "BR")
#' @param spreadsheet_path Full path to .csv, .xlsx, or .xls file 
#'        (required if species_input = "spreadsheet")
#' @param spreadsheet_sheet Name of Excel sheet/tab to read from 
#'        (ignored for CSV files)
#' @param species_column Name of column containing scientific names in the file
#' @param n_species Number of species to process (NULL = all species, useful for testing)
#' @param groups_column Column name containing group/taxon labels for filtering (optional)
#' @param selected_groups Character vector of group names to include (optional)
#' @param download_path Directory path where GBIF zip files will be saved 
#'        (default: current working directory)
#' @param output_filename Custom name for the downloaded zip file, without extension 
#'        (default: NULL, keeps GBIF's default key-based filename)
#' @param poll_interval Seconds between each status check while waiting 
#'        for GBIF to prepare the download (default: 10)
#' @return List containing: data (occurrence records), download_key, matching_results

download_gbif_occurrences <- function(species_input = NULL,
                                      taxon_group = NULL,
                                      rank = NULL,
                                      country_code = NULL,
                                      spreadsheet_path = NULL,
                                      spreadsheet_sheet = NULL,
                                      species_column = NULL,
                                      n_species = NULL,
                                      groups_column = NULL,
                                      selected_groups = NULL,
                                      download_path = NULL,
                                      output_filename = NULL,
                                      poll_interval = 10) {
  
  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (is.null(species_input) && is.null(taxon_group)) {
    stop("Either 'species_input' or 'taxon_group' must be provided.")
  }
  
  if (!is.null(species_input) && !is.null(taxon_group)) {
    warning("Both 'species_input' and 'taxon_group' provided. Using 'species_input'.")
    taxon_group <- NULL
  }
  
  if (!is.null(country_code) && !country_code %in% rgbif::enumeration_country()$iso2) {
    stop("Invalid country code: ", country_code)
  }
  
  if (is.null(download_path)) {
    download_path <- getwd()
  }
  
  if (!dir.exists(download_path)) {
    stop("Download path does not exist: ", download_path)
  }
  
  # ---------------------------------------------------------------------------
  # Process species input
  # ---------------------------------------------------------------------------
  if (!is.null(species_input)) {
    message("Processing species list...")
    species_list <- process_species_input(species_input, spreadsheet_path,
                                          spreadsheet_sheet, species_column,
                                          groups_column, selected_groups, n_species)
    download_type <- "species_list"
  } else {
    message("Processing taxonomic group: ", taxon_group)
    if (!is.null(rank)) message("Using rank: ", rank)
    species_list <- taxon_group
    download_type <- "taxon_group"
  }
  
  # ---------------------------------------------------------------------------
  # Match species to GBIF backbone taxonomy
  # ---------------------------------------------------------------------------
  message("Matching species names to GBIF backbone taxonomy...")
  gbif_keys_result <- get_gbif_taxonomic_keys(species_list, download_type, rank)
  
  if (length(gbif_keys_result$valid_keys) == 0) {
    stop("No valid taxonomic keys found. Check species names or taxonomic group.")
  }
  
  report_matching_results(gbif_keys_result, species_list)
  
  # ---------------------------------------------------------------------------
  # Build and submit download request
  # ---------------------------------------------------------------------------
  message("Submitting download request to GBIF...")
  
  if (download_type == "species_list") {
    if (!is.null(country_code)) {
      download_request <- rgbif::occ_download(
        rgbif::pred_in("taxonKey", gbif_keys_result$valid_keys),
        rgbif::pred("country", country_code),
        rgbif::pred("hasCoordinate", TRUE)
      )
    } else {
      download_request <- rgbif::occ_download(
        rgbif::pred_in("taxonKey", gbif_keys_result$valid_keys),
        rgbif::pred("hasCoordinate", TRUE)
      )
    }
  } else {
    if (!is.null(country_code)) {
      download_request <- rgbif::occ_download(
        rgbif::pred("taxonKey", gbif_keys_result$valid_keys[1]),
        rgbif::pred("country", country_code),
        rgbif::pred("hasCoordinate", TRUE)
      )
    } else {
      download_request <- rgbif::occ_download(
        rgbif::pred("taxonKey", gbif_keys_result$valid_keys[1]),
        rgbif::pred("hasCoordinate", TRUE)
      )
    }
  }
  
  download_key <- as.character(download_request)
  message("Download request created. GBIF key: ", download_key)
  
  # ---------------------------------------------------------------------------
  # Wait for GBIF to prepare the download (with live status updates)
  # ---------------------------------------------------------------------------
  message("\nWaiting for GBIF to prepare your download...")
  message(strrep("-", 50))
  
  start_time <- proc.time()["elapsed"]
  status <- "PREPARING"
  spinner <- c("|", "/", "-", "\\")
  tick <- 0
  
  while (!status %in% c("SUCCEEDED", "KILLED", "FAILED", "CANCELLED")) {
    
    Sys.sleep(poll_interval)
    tick <- tick + 1
    
    # Query current download status from GBIF
    meta <- tryCatch({
      rgbif::occ_download_meta(download_request)
    }, error = function(e) {
      message("Warning: could not fetch status - ", e$message)
      return(NULL)
    })
    
    if (!is.null(meta)) {
      status <- meta$status
    }
    
    # Calculate elapsed time
    elapsed <- round(proc.time()["elapsed"] - start_time)
    mins <- floor(elapsed / 60)
    secs <- elapsed %% 60
    elapsed_str <- sprintf("%02d:%02d", mins, secs)
    
    # Print live status line
    spin <- spinner[(tick %% length(spinner)) + 1]
    message(sprintf("  %s  Status: %-12s | Elapsed: %s | Attempt: %d",
                    spin, status, elapsed_str, tick))
  }
  
  message(strrep("-", 50))
  
  # Handle failed downloads
  if (status %in% c("KILLED", "FAILED", "CANCELLED")) {
    stop("GBIF download failed with status: ", status,
         "\nCheck your GBIF account at https://www.gbif.org/user/download for details.")
  }
  
  total_elapsed <- round(proc.time()["elapsed"] - start_time)
  message(sprintf("Download ready! Total wait time: %02d:%02d\n",
                  floor(total_elapsed / 60), total_elapsed %% 60))
  
  # ---------------------------------------------------------------------------
  # Download and (optionally) rename the zip file
  # ---------------------------------------------------------------------------
  message("Downloading zip file from GBIF servers...")
  gbif_zip <- rgbif::occ_download_get(download_request, path = download_path)
  
  # ---------------------------------------------------------------------------
  # Import occurrence data BEFORE renaming
  # ---------------------------------------------------------------------------
  message("Importing occurrence data...")
  gbif_data <- gbif_zip %>% rgbif::occ_download_import()
  
  # ---------------------------------------------------------------------------
  # Rename zip file AFTER import
  # ---------------------------------------------------------------------------
  if (!is.null(output_filename)) {
    
    # Ensure .zip extension
    if (!grepl("\\.zip$", output_filename)) {
      output_filename <- paste0(output_filename, ".zip")
    }
    
    # Construct original path from download key (GBIF always names it <key>.zip)
    original_path <- file.path(download_path, paste0(download_key, ".zip"))
    new_path      <- file.path(download_path, output_filename)
    
    if (file.exists(original_path)) {
      file.rename(original_path, new_path)
      message("Zip file renamed to: ", new_path)
    } else {
      warning("Could not rename zip file: original not found at ", original_path,
              "\nFile remains with GBIF default name.")
    }
  }
  
  report_download_results(gbif_data)
  
  return(list(
    data             = gbif_data,
    download_key     = download_key,
    matching_results = gbif_keys_result
  ))
}






# 2.2. Process Species Input from Various Sources ----------------------------
#' Handles different types of species input: direct character vectors, Excel 
#' spreadsheets (.xlsx, .xls), or CSV files (.csv).
#' Validates input sources and formats species lists for GBIF queries.
#'
#' @param species_input Either "spreadsheet" (works for .csv, .xlsx, and .xls) 
#'        or a character vector of species names
#' @param spreadsheet_path Path to the file (required if species_input = "spreadsheet")
#' @param spreadsheet_sheet Name of the sheet with species names 
#'        (only used for Excel files, ignored for CSV)
#' @param species_column Column name containing species names  
#' @param groups_column Column name containing group/taxon labels for filtering (optional)
#' @param selected_groups Character vector of group names to include (optional)
#' @param n_species Limit number of species to process (for testing)
#' @return Character vector of cleaned species names ready for GBIF query

process_species_input <- function(species_input,
                                  spreadsheet_path = NULL,
                                  spreadsheet_sheet = NULL,
                                  species_column = NULL,
                                  groups_column = NULL,              
                                  selected_groups = NULL,            
                                  n_species = NULL) {
  
  # ---------------------------------------------------------------------------
  # OPTION 1: Direct character vector input
  # ---------------------------------------------------------------------------
  if (is.character(species_input) && !identical(species_input, "spreadsheet")) {
    species_list <- species_input
    message("Species list provided directly as character vector.")
    
  # ---------------------------------------------------------------------------
  # OPTION 2: File input (.csv, .xlsx, or .xls)
  # ---------------------------------------------------------------------------
  } else if (identical(species_input, "spreadsheet") && !is.null(spreadsheet_path)) {
    
    # Validate species column is specified
    if (is.null(species_column)) {
      stop("'species_column' must be specified when reading from a file.")
    }
    
    # Read file (CSV or Excel) using helper function
    species_data <- read_species_file(spreadsheet_path, sheet = spreadsheet_sheet)
    
    # Validate species column exists
    if (!species_column %in% names(species_data)) {
      stop("Column '", species_column, "' not found in file. ",
           "Available columns: ", paste(names(species_data), collapse = ", "))
    }
    
    # Optional: filter by taxonomic group
    if (!is.null(groups_column) && !is.null(selected_groups)) {
      
      if (!groups_column %in% names(species_data)) {
        stop("Column '", groups_column, "' not found in file. ",
             "Available columns: ", paste(names(species_data), collapse = ", "))
      }
      
      available_groups <- unique(species_data[[groups_column]])
      available_groups <- available_groups[!is.na(available_groups)]
      message("Available groups in dataset: ", paste(available_groups, collapse = ", "))
      
      missing_groups <- setdiff(selected_groups, available_groups)
      if (length(missing_groups) > 0) {
        warning("The following selected groups were not found in the dataset: ", 
                paste(missing_groups, collapse = ", "))
      }
      
      valid_groups <- intersect(selected_groups, available_groups)
      if (length(valid_groups) == 0) {
        stop("None of the selected groups were found in the dataset.")
      }
      
      message("Filtering for groups: ", paste(valid_groups, collapse = ", "))
      species_data <- species_data %>%
        dplyr::filter(.data[[groups_column]] %in% valid_groups)
      
      group_summary <- species_data %>%
        dplyr::group_by(.data[[groups_column]]) %>%
        dplyr::summarise(n_species = dplyr::n(), .groups = "drop")
      
      message("Species count by selected groups:")
      for (i in 1:nrow(group_summary)) {
        message("  - ", group_summary[[groups_column]][i], ": ", 
                group_summary$n_species[i], " species")
      }
    }
    
    # Extract and clean species names
    species_list <- species_data[[species_column]]
    species_list <- species_list[!is.na(species_list) & species_list != ""]
    
    if (length(species_list) == 0) {
      stop("No valid species names found in the '", species_column, "' column.")
    }
    
  # ---------------------------------------------------------------------------
  # OPTION 3: Invalid input
  # ---------------------------------------------------------------------------
  } else {
    stop("Invalid 'species_input'. Provide either:\n",
         "  - A character vector of species names, or\n",
         "  - The keyword 'spreadsheet' along with 'spreadsheet_path' ",
         "(.csv, .xlsx, or .xls).")
  }
  
  # Optionally limit number of species for testing purposes
  if (!is.null(n_species) && n_species > 0) {
    if (n_species < length(species_list)) {
      message("Limiting to the first ", n_species, " species (out of ", 
              length(species_list), " total).")
      species_list <- head(species_list, n_species)
    }
  }
  
  message("Species list prepared: ", length(species_list), " species.")
  message("First few species: ", paste(head(species_list, 3), collapse = ", "))
  
  return(species_list)
}


# 2.3. Get GBIF Taxonomic Keys for Species Names ----------------------------
#' Queries GBIF backbone taxonomy to convert species names into taxonomic keys.
#' Handles both individual species lists and taxonomic groups with optional rank specification.
#'
#' @param species_list Character vector of species names
#' @param download_type Either "species_list" or "taxon_group"
#' @param rank Optional taxonomic rank for taxon_group (e.g., "class", "order", "family")
#' @return A list with:
#'   - valid_keys: vector of GBIF usageKeys
#'   - backbone_result: tibble with all matching info
#'   - unmatched_species: character vector of species not matched

get_gbif_taxonomic_keys <- function(species_list, download_type, rank = NULL) {
  
  if (download_type == "species_list") {
    message("Querying GBIF backbone for ", length(species_list), " species...")
    
    cat("Species to query:\n")
    for (i in seq_along(species_list)) {
      cat("  ", i, ": '", species_list[i], "'\n", sep = "")
    }
    
    backbone_result <- tryCatch({
      rgbif::name_backbone_checklist(species_list)
    }, error = function(e) {
      message("Error in name_backbone_checklist: ", e$message)
      return(NULL)
    })
    
    return(process_backbone_result(backbone_result, species_list))
    
  } else {
    message("Querying GBIF backbone for taxonomic group: ", species_list)
    if (!is.null(rank)) {
      message("Using rank specification: ", rank)
    }
    
    backbone_result <- NULL
    
    if (!is.null(rank)) {
      backbone_result <- tryCatch({
        message("Attempting lookup with rank '", rank, "'...")
        rgbif::name_backbone(species_list, rank = rank)
      }, error = function(e) {
        message("Error with rank specification: ", e$message)
        return(NULL)
      })
      
      if (!is.null(backbone_result) && !is.na(backbone_result$usageKey) && 
          backbone_result$matchType != "NONE") {
        message("Successfully matched '", species_list, "' with rank '", rank, "'")
      } else {
        message("Rank-based lookup failed. Trying without rank...")
        backbone_result <- NULL
      }
    }
    
    if (is.null(backbone_result)) {
      backbone_result <- tryCatch({
        message("Attempting lookup without rank specification...")
        rgbif::name_backbone(species_list)
      }, error = function(e) {
        message("Error in fallback lookup: ", e$message)
        return(NULL)
      })
    }
    
    if (is.null(backbone_result) || is.na(backbone_result$usageKey) || 
        backbone_result$matchType == "NONE") {
      return(list(
        valid_keys = integer(0),
        backbone_result = backbone_result,
        unmatched_species = species_list
      ))
    }
    
    return(list(
      valid_keys = backbone_result$usageKey,
      backbone_result = backbone_result,
      unmatched_species = character(0)
    ))
  }
}


# 2.4. Process Backbone Result for Species Lists ----------------------------
#' Process the backbone result from GBIF API with comprehensive error handling
#'
#' @param backbone_result Data frame from GBIF API
#' @param species_list Original species list
#' @return List with processed results

process_backbone_result <- function(backbone_result, species_list) {
  
  if (is.null(backbone_result)) {
    message("backbone_result is NULL")
    return(list(
      valid_keys = integer(0),
      backbone_result = NULL,
      unmatched_species = species_list
    ))
  }
  
  message("Debug - backbone_result structure:")
  message("  Class: ", class(backbone_result))
  message("  Dimensions: ", paste(dim(backbone_result), collapse = " x "))
  message("  Column names: ", paste(names(backbone_result), collapse = ", "))
  
  if (nrow(backbone_result) > 0) {
    message("First few rows of backbone_result:")
    print(head(backbone_result, 3))
  }
  
  required_cols <- c("matchType", "usageKey")
  missing_cols <- required_cols[!required_cols %in% names(backbone_result)]
  
  if (length(missing_cols) > 0) {
    warning("Missing required columns: ", paste(missing_cols, collapse = ", "))
    
    if ("speciesKey" %in% names(backbone_result)) {
      message("Using 'speciesKey' instead of 'usageKey'")
      backbone_result$usageKey <- backbone_result$speciesKey
    } else if ("acceptedUsageKey" %in% names(backbone_result)) {
      message("Using 'acceptedUsageKey' instead of 'usageKey'")
      backbone_result$usageKey <- backbone_result$acceptedUsageKey
    } else {
      message("No suitable key column found")
      return(list(
        valid_keys = integer(0),
        backbone_result = backbone_result,
        unmatched_species = species_list
      ))
    }
  }
  
  if (!"matchType" %in% names(backbone_result)) {
    message("matchType column not found. Checking for alternatives...")
    if ("status" %in% names(backbone_result)) {
      message("Using 'status' column as matchType proxy")
      backbone_result$matchType <- backbone_result$status
    } else {
      message("Creating default matchType based on usageKey presence")
      backbone_result$matchType <- ifelse(is.na(backbone_result$usageKey), "NONE", "EXACT")
    }
  }
  
  tryCatch({
    has_usage_key <- !is.na(backbone_result$usageKey) & backbone_result$usageKey != ""
    has_valid_match <- backbone_result$matchType != "NONE" & !is.na(backbone_result$matchType)
    valid_rows <- has_usage_key & has_valid_match
    
    message("Debug - filtering:")
    message("  Total rows: ", nrow(backbone_result))
    message("  Rows with usageKey: ", sum(has_usage_key, na.rm = TRUE))
    message("  Rows with valid matchType: ", sum(has_valid_match, na.rm = TRUE))
    message("  Valid rows: ", sum(valid_rows, na.rm = TRUE))
    
    if (sum(valid_rows, na.rm = TRUE) == 0) {
      message("No valid matches found")
      return(list(
        valid_keys = integer(0),
        backbone_result = backbone_result,
        unmatched_species = species_list
      ))
    }
    
    valid_matches <- backbone_result[valid_rows, ]
    
    species_only <- valid_matches %>%
      dplyr::filter(
        rank == "SPECIES",
        status == "ACCEPTED"
      )
    
    message("Debug - species filtering:")
    message("  Valid matches before filter: ", nrow(valid_matches))
    message("  Species-only matches after filter: ", nrow(species_only))
    
    if (nrow(species_only) > 0 && "speciesKey" %in% names(species_only)) {
      valid_keys <- species_only$speciesKey
      message("  Using speciesKey for taxonomic keys")
    } else {
      valid_keys <- species_only$usageKey
      message("  Fallback to usageKey (speciesKey not available)")
    }
    
    valid_keys <- valid_keys[!is.na(valid_keys)]
    
    unmatched_species <- character(0)
    if (nrow(species_only) > 0 && "verbatim_name" %in% names(species_only)) {
      matched_names <- species_only$verbatim_name
      unmatched_species <- setdiff(species_list, matched_names)
    } else if (nrow(species_only) > 0 && "scientificName" %in% names(species_only)) {
      matched_names <- species_only$scientificName
      unmatched_species <- setdiff(species_list, matched_names)
    } else {
      n_matched <- length(valid_keys)
      n_total <- length(species_list)
      if (n_matched < n_total) {
        unmatched_species <- species_list[(n_matched + 1):n_total]
      }
    }
    
    return(list(
      valid_keys = valid_keys,
      backbone_result = backbone_result,
      unmatched_species = unmatched_species
    ))
    
  }, error = function(e) {
    message("Error in filtering valid matches: ", e$message)
    message("Returning empty result")
    return(list(
      valid_keys = integer(0),
      backbone_result = backbone_result,
      unmatched_species = species_list
    ))
  })
}


# 2.5. Report Taxonomic Matching Results ------------------------------------
#' Displays summary statistics about how many species names were successfully
#' matched to GBIF backbone taxonomy vs. those that couldn't be resolved.
#'
#' @param gbif_keys_result Results from get_gbif_taxonomic_keys function
#' @param original_input Original species list or taxonomic group provided by user

report_matching_results <- function(gbif_keys_result, original_input) {
  
  total_input <- length(original_input)
  valid_matches <- length(gbif_keys_result$valid_keys)
  unmatched <- length(gbif_keys_result$unmatched_species)
  
  message("Taxonomic matching results:")
  message("  - Input species/groups: ", total_input)
  message("  - Successfully matched: ", valid_matches)
  message("  - Unmatched: ", unmatched)
  
  if (valid_matches > 0) {
    message("  - Valid GBIF keys found: ", 
            paste(head(gbif_keys_result$valid_keys, 5), collapse = ", "))
  }
  
  if (unmatched > 0 && length(gbif_keys_result$unmatched_species) > 0) {
    warning("The following species were not matched in GBIF backbone:")
    for (sp in head(gbif_keys_result$unmatched_species, 10)) {
      message("    - ", sp)
    }
    if (unmatched > 10) {
      message("    ... and ", unmatched - 10, " more")
    }
  }
}


# 2.6. Report Download Results Summary --------------------------------------
#' Displays summary statistics about the downloaded occurrence data including
#' total records, unique species, and top species by occurrence count.
#'
#' @param gbif_data Downloaded occurrence dataset from GBIF

report_download_results <- function(gbif_data) {
  
  n_records <- nrow(gbif_data)
  n_species <- length(unique(gbif_data$species))
  
  message("Download completed successfully!")
  message("  - Total records: ", format(n_records, big.mark = ","))
  message("  - Unique species: ", n_species)
  
  if (n_species > 0) {
    species_counts <- gbif_data %>%
      dplyr::count(species, sort = TRUE)
    
    message("  - Top species by occurrence count:")
    for (i in 1:nrow(species_counts)) {
      message("    ", i, ". ", species_counts$species[i], ": ", 
              format(species_counts$n[i], big.mark = ","))
    }
  }
}


# 2.7. DEBUG FUNCTION - Test name resolution without downloading -------------
#' Standalone utility to troubleshoot species name matching against the GBIF 
#' backbone taxonomy without triggering a full download.
#' Supports direct character vectors, .csv, .xlsx, and .xls files.
#'
#' When to use:
#' ✅ Diagnostic tool – check why species names are not matched by GBIF  
#' ✅ Quick testing – validate a list of names before running a full download  
#' ✅ Understand API responses – see exactly what GBIF returns for each query  
#' ✅ Development/debugging – useful while modifying the main pipeline  
#'
#' @param species_input Either "spreadsheet" (works for .csv, .xlsx, and .xls) 
#'        or a character vector of species names
#' @param spreadsheet_path Path to the file (required if species_input = "spreadsheet")
#' @param spreadsheet_sheet Name of the sheet with species names 
#'        (only used for Excel files, ignored for CSV)
#' @param species_column Column name containing species names
#' @param groups_column Column name containing group/taxon labels for filtering (optional)
#' @param selected_groups Character vector of group names to include (optional)
#' @param n_species Number of species to test (NULL = all; keep small for debugging)
#' @return List with valid_keys, backbone_result, and unmatched_species

debug_name_matching <- function(species_input,
                                spreadsheet_path = NULL,
                                spreadsheet_sheet = NULL,
                                species_column = NULL,
                                groups_column = NULL,
                                selected_groups = NULL,
                                n_species = NULL) {
  
  message("=== DEBUGGING NAME MATCHING ===")
  
  # Use process_species_input which now handles CSV, xlsx, and xls via read_species_file
  species_list <- process_species_input(
    species_input     = species_input,
    spreadsheet_path  = spreadsheet_path,
    spreadsheet_sheet = spreadsheet_sheet,
    species_column    = species_column,
    groups_column     = groups_column,
    selected_groups   = selected_groups,
    n_species         = n_species
  )
  
  # ---------------------------------------------------------------------------
  # Test individual species lookups (first 10 only)
  # ---------------------------------------------------------------------------
  message("\n=== TESTING INDIVIDUAL LOOKUPS ===")
  for (i in seq_along(head(species_list, 10))) {
    sp <- species_list[i]
    message("\nTesting species ", i, ": '", sp, "'")
    
    individual_result <- tryCatch({
      rgbif::name_backbone(sp)
    }, error = function(e) {
      message("  Error: ", e$message)
      return(NULL)
    })
    
    if (!is.null(individual_result)) {
      message("  Match type: ", individual_result$matchType)
      message("  Usage key: ", individual_result$usageKey)
      message("  Scientific name: ", individual_result$scientificName)
    }
  }
  
  # ---------------------------------------------------------------------------
  # Test batch lookup
  # ---------------------------------------------------------------------------
  message("\n=== TESTING BATCH LOOKUP ===")
  batch_result <- get_gbif_taxonomic_keys(species_list, "species_list")
  
  return(batch_result)
}
