# wmhelpers


**wmhelpers** is an internal R package developed by **WildMon (WM)** to support and standardize ecological analysis workflows across projects and collaborators.

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
>
> ```r
> ?read_issue_sites
> ```

---

## 📖 Documentation philosophy

* The **README** provides a high-level overview of the package and its scope
* **Each function is fully documented** using Roxygen2 and accessible via `?function_name`
* Detailed examples, parameters, and edge cases live in the function-level documentation

This approach keeps the README concise while ensuring technical details remain discoverable and up to date.

---

## 👥 Contributors

* **Jose W. Ribeiro Jr.** ("Xuleta")
* **Thiago Bicudo**

Additional contributors should be added via GitHub pull requests.

---

## ⚙️ Development status

`wmhelpers` is under active development. New helpers will be added incrementally as workflows evolve.

Breaking changes are avoided whenever possible and documented when unavoidable.

---

## 📄 License

This project is released under the MIT License.

---

## 📬 Contact

For questions, suggestions, or bug reports, please open an issue on GitHub or contact the WildMon team.
