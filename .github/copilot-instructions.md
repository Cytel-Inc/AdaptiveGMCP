# AI Coding Agent Instructions — AdaptiveGMCP

These instructions make AI agents immediately productive in this R package. They document actual patterns and workflows observed in the repo.

## Overview
- Purpose: Graph-based, adaptive multi-arm, multi-endpoint clinical trial methods with simulation and analysis.
- Key exports: `simMAMSMEP`, `adaptGMCP_CER`, `adaptGMCP_PC`, `plotGraph`, `genPowerTablePlots`, `AdaptGMCPSimApp` (Shiny).
- See exports in [NAMESPACE](../NAMESPACE) and package metadata in [DESCRIPTION](../DESCRIPTION).

## Architecture
- Simulation flow: `simMAMSMEP` builds a `gmcpSimObj` and drives per-look simulation and selection; see [R/MAMSMEP_SIMULATION_MAIN.R](../R/MAMSMEP_SIMULATION_MAIN.R) and wrapper in [R/simMAMSMEP_Wrapper.R](../R/simMAMSMEP_Wrapper.R).
- Analysis flows:
  - CER method: [R/cerAdaptGMCP_Analysis.R](../R/cerAdaptGMCP_Analysis.R) — stage 1/2 analysis with conditional error computation and optional adaptation.
  - P‑value combination: [R/pValueAdaptGMCP_Analysis.R](../R/pValueAdaptGMCP_Analysis.R) — inverse normal combination with rpact boundaries.
- Graph visualization: [R/graphPlot.R](../R/graphPlot.R) renders visNetwork graphs of hypotheses, weights, and transitions.
- Shiny app: [inst/shinyApps/AdaptGMCPSimApp.R](../inst/shinyApps/AdaptGMCPSimApp.R) with modules in the same folder (e.g., `tabularModule.R`, `correlationMatrixModule.R`).

## Data & Hypothesis Conventions
- Hypothesis ordering: iterate endpoints within each treatment before the next endpoint (e.g., for 2 trts × 3 eps: H1=EP1/T1, H2=EP1/T2, H3=EP2/T1, H4=EP2/T2, H5=EP3/T1, H6=EP3/T2). See UI note in [inst/shinyApps/AdaptGMCPSimApp.R](../inst/shinyApps/AdaptGMCPSimApp.R) and `WI`/`G` usage across functions.
- Weights and transitions: `WI` (initial weights) and `G` (transition matrix) drive intersection weights via `genWeights()`; conventions are consistent across simulation and analyses.
- Info fraction and allocation: `info_frac` is cumulative per-look; `simMAMSMEP` normalizes `Arms.alloc.ratio` so control equals 1.
- Two-arm constraint: For 2 arms, parametric tests are auto-downgraded (Bonferroni / Non‑Parametric); see logic in [R/MAMSMEP_SIMULATION_MAIN.R](../R/MAMSMEP_SIMULATION_MAIN.R).

## Build & Dev Workflow
- Helper script: [rebuild_package.R](../rebuild_package.R) exposes common tasks:
  - Quick dev loop:
    ```r
    source('rebuild_package.R'); quick_reinstall(); load_pkg()
    ```
  - Full rebuild:
    ```r
    source('rebuild_package.R'); full_rebuild()
    ```
  - Docs & checks:
    ```r
    source('rebuild_package.R'); document_pkg(); check_pkg(); test_pkg()
    ```
  - Dependencies:
    ```r
    source('rebuild_package.R'); install_deps()
    ```
- Roxygen: Functions are documented inline; run `document_pkg()` to refresh man pages and `NAMESPACE`.

## Tests
- Test harness: [tests/testthat.R](../tests/testthat.R) with cases in [tests/testthat](../tests/testthat).
- Run tests:
  ```r
  source('rebuild_package.R'); test_pkg()
  # or
  devtools::test()
  ```
- Snapshot files live under [tests/testthat/_snaps](../tests/testthat/_snaps).

## Linting & Quality
- Lintr script: [internalData/LintThisPkg.R](../internalData/LintThisPkg.R) shows typical usage and exclusions.
- Run lint:
  ```r
  lintr::lint_dir('R')
  # or source the helper and use targeted linters
  ```

## Shiny App
- Entry: `AdaptGMCPSimApp` (exported; UI/server in [inst/shinyApps/AdaptGMCPSimApp.R](../inst/shinyApps/AdaptGMCPSimApp.R)).
- Modules: `tabularModule.R`, `correlationMatrixModule.R`, `transitionMatrixModule.R`, `IAModule.R`, `helper_inputDataCsv.R` under the same folder.
- Typical launch:
  ```r
  library(AdaptGMCP); AdaptGMCPSimApp()
  ```

## External Dependencies
- rpact: group sequential design boundaries and alpha spending (used to derive `stageLevels` and `alphaSpent`).
- Core runtime: `data.table`, `dplyr`, `mvtnorm`, `Matrix`, `matrixcalc`, `parallel`, `stringr`.
- Visualization/UI: `visNetwork`, `ggplot2`, `gridExtra`, `shiny`, `rhandsontable`, `htmltools`, `shinycssloaders`.
- See [DESCRIPTION](../DESCRIPTION) for full `Imports`/`Suggests`.

## Agent Tips
- Preserve exported signatures and documented defaults; they are referenced by examples in [internalData](../internalData).
- Keep hypothesis ordering consistent when generating or consuming `WI`, `G`, and tables.
- For boundaries, prefer rpact-derived `userAlphaSpending` unless a domain change requires custom spending.
- Use `renv::restore()` before running tests locally to match lockfile; helper functions in `rebuild_package.R` wrap common renv actions.

## Style Guide
- Follow the Google R Style Guide for all new code and significant changes: https://google.github.io/styleguide/Rguide.html (Cytel preference).
- Do not refactor legacy code solely to satisfy style; limit changes to scope of the task.
- Basics to keep consistent: two-space indentation, use `<-` for assignment, `TRUE`/`FALSE` (not `T`/`F`), spaces around operators and after commas, and avoid `attach()`.
- Naming conventions (Google guide):
  - Variables: lowercase with dots (e.g., `variable.name`)
  - Functions: UpperCamelCase for exported functions, lowerCamelCase for internal helpers
  - Constants: `kConstantName`
  - Files: lowercase with underscores or hyphens
- Use roxygen2 for documentation; after changes run `document_pkg()` to refresh man pages and [NAMESPACE](../NAMESPACE).
- Lint with `lintr::lint_dir('R')`; see [internalData/LintThisPkg.R](../internalData/LintThisPkg.R) for targeted linters and exclusions.
