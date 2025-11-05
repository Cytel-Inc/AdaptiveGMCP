# LintThisPkg.R
# Code in this file is used for linting the code in this package.
library(lintr)

### USING .lintr FILE ###########################
  # Lint everything first
  all_lints <- lintr::lint_dir("R")

  # Filter out certain files
  excluded <- c("AdaptGMCPSimulationShinyApp.R",
                "graphicalMCP_pow_calc_wrapper.R",
                "gMCP_calcPower_wrapper.R")
  lints <- all_lints[!vapply(all_lints, \(x) x$filename %in% excluded, logical(1))]

  # Inspect or summarize remaining results
  unique(vapply(lints, `[[`, character(1), "filename"))

  lints
################################################

### RUNNING SPECIFIC LINTERS ###################
# This is bypassing the .lintr file.
  lintr::lint_dir(path = "R", linters = lintr::cyclocomp_linter(15))

  lintr::lint_dir(path = "R", linters = lintr::fixed_regex_linter(15))
  ################################################
