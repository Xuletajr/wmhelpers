# Development Guidelines and Best Practices

Thank you for contributing to this package.

This project follows strict development standards to ensure consistency, reproducibility, and long-term maintainability. All contributors must adhere to the guidelines below.

---

## 1. General Development Philosophy

This package is designed to be:

- Reproducible  
- Scientifically robust  
- Scalable across biodiversity monitoring projects  
- Maintainable by multiple contributors  

Consistency, clarity, and documentation are mandatory — not optional.

---

## 2. Function Naming Convention

All function names must follow a **verb + object** structure.

### General Pattern

- `read_<data_type>_<object>()`
- `build_<object>()`
- `create_<object>()`
- `validate_<object>()`
- `export_<object>()`

### Examples

- `read_camtrap_detections()`
- `read_audio_detections()`
- `build_y_array()`
- `validate_deployment_table()`

### Rules

- Use snake_case.
- Use lowercase.
- Start with a clear action verb.
- Avoid vague names like `process_data()` or `run_model()`.
- Do not duplicate existing functionality under a new name.

Before creating a new function, verify that similar logic does not already exist in the package.

---

## 3. Internal vs Exported Functions

### Exported functions

- Must include `@export` in roxygen2.
- Must be fully documented.
- Must follow naming standards.

### Internal helper functions

- Should start with a leading dot:

```r
.build_detection_matrix()
```

- Should NOT include `@export`.
- Are intended for internal use only.

---

## 4. Mandatory roxygen2 Documentation

All exported functions must include complete roxygen2 headers.

Each function must contain:

- `@param` for every argument  
- `@return`  
- `@details` (when appropriate)  
- `@export`  

### Required Structure

```r
#' Title
#'
#' Description of what the function does.
#'
#' @param argument_name Description of the argument.
#' @return Description of the returned object.
#' @details Additional explanation if needed.
#'
#' @export
function_name <- function(argument_name) {

  # function body

}
```

---

## 5. Required Workflow Before Commit

After:

- Adding a new function  
- Modifying function arguments  
- Editing documentation  
- Changing exported functions  

You must run:

```r
devtools::document()
```

Then run:

```r
devtools::check()
```

Ensure no new errors or warnings are introduced before committing.

---

## 6. Pre-Commit Checklist

Before committing and pushing changes, ensure that:

- [ ] The function follows the established naming conventions (verb + object).
- [ ] No duplicated logic exists elsewhere in the package.
- [ ] All parameters are fully documented using roxygen2.
- [ ] The function includes `@export` when appropriate.
- [ ] `devtools::document()` has been executed.
- [ ] `devtools::check()` runs without introducing new errors or warnings.

Only after completing this checklist should changes be committed and pushed.

---

## 7. Pull Requests

When submitting a pull request:

- Provide a clear and concise description of what was changed.
- Explain why the change is necessary.
- Keep pull requests focused and minimal (avoid mixing unrelated changes).
- Ensure the branch is up to date with `main`.
- Confirm that all development and documentation standards described in this document have been followed.

Pull requests that do not comply with these guidelines may be requested for revision before approval.
