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

### 📁 Core data I/O and detection

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

### 📊 Visualization

| Function | Description |
|---|---|
| `plot_bar_species()` | Generate standardized bar plots for species detection frequency and richness summaries |
	


---

### 🌐 GBIF occurrence download

Functions to query, match, and download species occurrence data from the [GBIF](https://www.gbif.org) database.

| Function | Description |
|---|---|


---

### 🦉 MSOM — General inference and diagnostics (`msom_general.R`)

Functions for MCMC diagnostics, species-level parameter summaries, and occupancy probability extraction from fitted Multi-Species Occupancy Models (output of package spOccupancy).

#### Diagnostics

| Function | Description |
|---|---|


#### Parameter summaries

| Function | Description |
|---|---|


#### Detection prediction over time

| Function | Description |
|---|---|


#### Occupancy prediction pipeline

| Function | Description |
|---|---|


#### Site-level occupancy (psi)

| Function | Description |
|---|---|

---

### 📈 MSOM — Covariate effects and marginal plots (`msom_cov_effects.R`)

Functions for summarising and visualising species-level marginal effects of environmental covariates.

| Function | Description |
|---|---|

---

### 🗺️ MSOM — Spatial prediction (`msom_predict.R`)

Functions for preparing covariates, reprojecting coordinates, and running large-scale spatial occupancy predictions on H3 hexagonal grids.

| Function | Description |
|---|---|

---

### 🗺️ MSOM — Prediction maps and rasters (`msom_prediction_maps.R`)

End-to-end functions for loading predictions, building spatial objects, and exporting publication-quality maps and GeoTIFFs.

| Function | Description |
|---|---|


---

### 📋 Report Appendix

| Function | Description |
|---|---|

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
