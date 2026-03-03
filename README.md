# wmhelpers
![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)
![R version](https://img.shields.io/badge/R-%E2%89%A5%204.2-blue.svg)
![GitHub last commit](https://img.shields.io/github/last-commit/WM/wmhelpers)
![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)
<p align="left">
  <img src="inst/figures/Wildmon_Horizontal_logo_white_blue_full.png"
       alt="WildMon logo"
       width="220"/>
</p>

**wmhelpers** is an internal R package developed by **WildMon (WM)** Quantitative Ecology Team to support and standardize ecological analysis workflows across projects and collaborators.

The package provides a growing collection of **helper functions** to manipulate, prepare, quality-check, and automate data processing steps commonly required in ecological analyses using data derived from multiple monitoring technologies.

---

## 🎯 Purpose

The main goals of `wmhelpers` are to:

* Centralize reusable helper functions used across WildMon R workflows
* Reduce code duplication in analysis scripts
* Improve reproducibility and consistency across projects
* Minimize errors caused by diverging local implementations
* Facilitate collaborative development and maintenance of shared logic

This package is **not** intended to perform full analyses, but rather to support them.

---

## 🌱 Supported data sources

Functions in `wmhelpers` are designed to assist workflows involving ecological data from, but not limited to:

* Passive Acoustic Monitoring (PAM)
* Camera traps
* Environmental DNA (eDNA)
* Citizen science platforms
* GBIF biodiversity occurrence databases
* Multi-Species Occupancy Models (MSOM)
* Project tracking and metadata tables

---

## 📦 Installation

`wmhelpers` is an internal package hosted on GitHub.

```r
remotes::install_github("WM/wmhelpers")
```

Load the package:

```r
library(wmhelpers)
```

---

## 🧰 Available functions

Functions are organized by module. All exported functions are documented and accessible via R help (e.g. `?function_name`).

---

### 📁 Core data I/O and detection helpers

| Function | Description |
|---|---|
| `build_detection_array()` | Construct a 3D detection array (species × sites × dates) where 1 = detection, 0 = sampled with no detection, and NA = no sampling effort |
| `coalesce_suffix_pairs()` | Merge duplicated columns created by multi-method joins using configurable suffix priority |
| `read_audio_detections()` | Read and clean confirmed acoustic detections within deployment windows |
| `read_camtrap_detections()` | Read and clean camera trap detections from all `images*.csv` files |
| `read_deployment_file()` | Read and standardize deployment metadata (start/end dates per site) |
| `read_issue_sites()` | Read and standardize site exclusion flags from a Project Tracking Google Sheets document |
| `read_species_list()` | Read, clean, and standardize species metadata from a Project Tracking sheet |
| `validate_detection_array()` | Validate consistency between detection arrays and raw multi-method detection data |

---

### 📊 Visualization helpers

| Function | Description |
|---|---|
| `plot_bar_species()` | Generate standardized bar plots for species detection frequency and richness summaries |
| `bar_chart()` | Automated bar chart with configurable IUCN status labels and Species of Greatest Conservation Need (SGCN) flags |

---

### 🌐 GBIF occurrence download

Functions to query, match, and download species occurrence data from the [GBIF](https://www.gbif.org) database.

| Function | Description |
|---|---|
| `download_gbif_occurrences()` | Main wrapper: download species occurrence records from GBIF, supporting species lists (from vectors, CSV, or Excel) and broad taxonomic groups, with live status polling |
| `process_species_input()` | Parse species input from direct character vectors, `.csv`, `.xlsx`, or `.xls` files, with optional group filtering |
| `get_gbif_taxonomic_keys()` | Match species names or taxonomic groups to the GBIF backbone taxonomy and return valid usage keys |
| `process_backbone_result()` | Process and validate the raw backbone API result, with fallback logic for non-standard column structures |
| `report_matching_results()` | Print a summary of how many species were successfully matched vs. unresolved |
| `report_download_results()` | Print a summary of the downloaded dataset (total records, unique species, top occurrences) |
| `debug_name_matching()` | Diagnostic utility: test species name resolution without triggering a full download; supports vectors, CSV, and Excel |

---

### 🦉 MSOM — General inference and diagnostics (`msom_general.R`)

Functions for MCMC diagnostics, species-level parameter summaries, and occupancy probability extraction from fitted Multi-Species Occupancy Models.

#### Diagnostics

| Function | Description |
|---|---|
| `check_rhat()` | Evaluate R-hat (Gelman-Rubin) convergence per parameter and return a table of problematic ones |
| `plot_rhat_hist()` | Faceted histogram of R-hat values for beta (occupancy) and alpha (detection) parameters |
| `check_ess()` | Identify parameters with Effective Sample Size (ESS) below a minimum threshold |
| `plot_ess_hist()` | Faceted histogram of ESS values for beta and alpha parameters |

#### Parameter summaries

| Function | Description |
|---|---|
| `extract_prob()` | Back-transform MCMC samples to the probability scale and summarise occupancy or detection per species and class/year |
| `build_param_table()` | Build a tidy species-level summary table (mean, SD, 95% CRI, effect direction) from beta or alpha samples |
| `get_significant_vars()` | Identify community-level variables whose 95% CRI entirely excludes zero |
| `plot_occ_vs_det()` | Scatter plot of mean occupancy vs. mean detection probability with linear regression annotation |

#### Detection prediction over time

| Function | Description |
|---|---|
| `predict_detection_julian()` | Predict community-level detection probability across Julian days using linear + quadratic effects |
| `plot_detection_julian()` | Ribbon + line plot of community-level detection probability over calendar survey dates |

#### Occupancy prediction pipeline

| Function | Description |
|---|---|
| `prepare_occ_sequences()` | Generate raw and scaled prediction sequences for each continuous covariate |
| `predict_occupancy()` | Fill a 3D prediction array `[n_pred × nsamp × combinations]` from community-level MCMC beta samples |
| `run_occ_predictions()` | Master wrapper that sequentially calls `prepare_occ_sequences()`, `predict_occupancy()`, and `summarise_occ_predictions()` |
| `plot_occ_predictions()` | Produce a named list of ggplot objects (one per variable × year × class combination) |
| `export_occ_plots()` | Save individual PNGs and an optional multiplot for selected variables |

#### Site-level occupancy (psi)

| Function | Description |
|---|---|
| `extract_psi_stats()` | Compute mean, SD, and 95% CRI from the posterior `psi.samples` array per species × site |
| `export_psi_tables()` | Write one CSV per statistic (mean, SD, CRI 0.025, CRI 0.975) in sites × species format |

---

### 📈 MSOM — Covariate effects and marginal plots (`msom_cov_effects.R`)

Functions for summarising and visualising species-level marginal effects of environmental covariates.

| Function | Description |
|---|---|
| `get_xlabel()` | Return a human-readable x-axis label for a covariate name, with support for custom overrides |
| `clean_covariate_name()` | Sanitise a covariate string (remove parentheses, carets, spaces) for safe use in file paths |
| `get_cov_info()` | Compute the CRI range and effect-direction count for a given covariate, grouped by a trait variable |
| `export_cov_info()` | Iterate over a covariate vector and export all summaries to a plain-text `.txt` file |
| `plot_covariate()` | Build, save, and return dot-and-errorbar plots of species-level marginal effects, with optional faceting and automatic bird-dataset splitting |
| `move_output_figures()` | Organise output PNGs into `1.species/` and `2.community/` subfolders with optional taxon sub-subfolders |

---

### 🗺️ MSOM — Spatial prediction (`msom_predict.R`)

Functions for preparing covariates, reprojecting coordinates, and running large-scale spatial occupancy predictions on H3 hexagonal grids.

| Function | Description |
|---|---|
| `load_hexgrid_variables()` | Auto-detect and load a hexgrid CSV, lowercase columns, remove NA rows, rename, and apply optional column transformations |
| `prepare_hex_covariates()` | Select model covariates, standardise, and build a design matrix with intercept ready for `predict()` |
| `project_coords()` | Reproject lon/lat from WGS84 to a user-specified metric CRS and return a tibble of projected coordinates |
| `predict_occupancy_chunks()` | Run chunked spatial prediction over a hexgrid, returning 3D arrays of species occupancy and community richness quantiles |
| `plot_species_occupancy_map()` | Quick-check choropleth map of posterior median occupancy for a single species on the H3 hexgrid |

---

### 🗺️ MSOM — Prediction maps and rasters (`msom_prediction_maps.R`)

End-to-end functions for loading predictions, building spatial objects, and exporting publication-quality maps and GeoTIFFs.

| Function | Description |
|---|---|
| `load_msom_data()` | Load the MSOM workspace RData file and return filtered species codes, scientific names, and common names |
| `create_msom_dirs()` | Create the output directory structure for maps and rasters |
| `import_msom_predictions()` | Read the RDS files containing psi and richness quantile arrays produced by the prediction pipeline |
| `load_detection_table()` | Read the raw species-by-site detection CSV and subset to modelled species |
| `load_site_shapefile()` | Read the site point shapefile, standardise names (with optional corrections), and filter to surveyed sites |
| `build_hex_grid()` | Convert H3 cell addresses from the prediction array into an `sf` hexagonal polygon grid |
| `build_utm_points()` | Project site locations, extract UTM coordinates, and join with the detection table |
| `plot_predictions()` | Generate and save a species occupancy map (hexgrid fill + detection points) as a PNG |
| `plot_richness_map()` | Generate and save a community richness map (posterior median) as a PNG |
| `export_species_rasters()` | Rasterise per-species median occupancy predictions and write individual GeoTIFF files |

---

### 📋 Reporting helpers

| Function | Description |
|---|---|
| `update.appendix()` | Join occupancy and detection probability estimates to an existing appendix data frame, supporting both intercept-only and multi-class/year model outputs |

---

> 📌 Detailed documentation for each function is available via R help (e.g. `?plot_covariate`, `?download_gbif_occurrences`).

---

## Contributing

Please read the [CONTRIBUTING.md](CONTRIBUTING.md) file before submitting changes.

All contributions must follow the development standards and documentation guidelines described there.

---

## 👥 Contributors

* **Jose W. Ribeiro Jr.** ("Xuleta")
* **Thiago Bicudo**

---

## ⚙️ Development status

`wmhelpers` is under active development. New helpers will be added incrementally as workflows evolve.

Breaking changes are avoided whenever possible and documented when unavoidable.

---

## 📬 Contact

For questions, suggestions, or bug reports, please open an issue on GitHub or contact the WildMon Quantitative Ecology team.
