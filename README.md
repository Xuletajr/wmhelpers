# wmhelpers

![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)
![R version](https://img.shields.io/badge/R-%E2%89%A5%204.2-blue.svg)
![GitHub last commit](https://img.shields.io/github/last-commit/WM/wmhelpers)
![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)


<p align="left">
  <img src="inst/figures/Wildmon-Horizontal-logo-white-blue-full.png"
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

The table below lists the currently available exported helpers.

| Function             | Description                                                                              |
| -------------------- | ---------------------------------------------------------------------------------------- |
| `read_issue_sites()` | Read and standardize site exclusion flags from a Project Tracking Google Sheets document |

> 📌 Detailed documentation for each function is available via R help:

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
