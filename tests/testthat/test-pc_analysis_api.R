############
# Test1
testthat::test_that("Test 1: PC analysis API scaffolds", {
  testthat::skip_if_not(
    isTRUE(getOption("AdaptGMCP.run_pc_analysis_api_tests", FALSE)),
    "Scaffold only: set options(AdaptGMCP.run_pc_analysis_api_tests = TRUE) and fill expected values"
  )

  wi <- c(1 / 2, 1 / 2, 0, 0)
  g <- matrix(c(
    0, 1 / 2, 1 / 2, 0,
    1 / 2, 0, 0, 1 / 2,
    0, 1, 0, 0,
    1, 0, 0, 0
  ), byrow = TRUE, nrow = 4)

  t <- c(0.5, 0.7, 1)
  alp <- 0.025

  corr <- matrix(c(
    1, 0.5, 0.5, NA,
    0.5, 1, NA, 0.5,
    0.5, NA, 1, 0.5,
    NA, 0.5, 0.5, 1
  ), byrow = TRUE, nrow = 4)

  tt <- "Partly-Parametric"
  des <- "asOF"

  state <- SetupAnalysis_PC(
    WI = wi,
    G = g,
    test.type = tt,
    alpha = alp,
    info_frac = t,
    typeOfDesign = des,
    Correlation = corr,
    plotGraphs = FALSE
  )

  testthat::expect_s3_class(state, "PCAnalysisState")
  testthat::expect_equal(state$completed_looks, 0L)

  # Look 1
  state <- AnalyzeLook_PC(
    state,
    p_raw = c(H1 = 0.01, H2 = 0.20, H3 = 0.15, H4 = 0.30),
    plotGraphs = FALSE
  )

  testthat::expect_equal(state$completed_looks, 1L)

  exp_out <- readRDS("t1.l1.AdjPValues.rds")
  expected_cols_look1 <- colnames(exp_out)
  testthat::expect_true(all(expected_cols_look1 %in% colnames(state$mcpObj$AdjPValues)))
  testthat::expect_equal(state$mcpObj$AdjPValues, exp_out)

  exp_out <- readRDS("t1.l1.RejFlagCurr.rds")
  testthat::expect_equal(state$mcpObj$rej_flag_Curr, exp_out)

  # Look 2 with selection
  state <- AnalyzeLook_PC(
    state,
    p_raw = c(H1 = 0.02, H2 = 0.10, H4 = 0.40),
    selection = c("H1", "H2", "H4"),
    plotGraphs = FALSE
  )

  testthat::expect_equal(state$completed_looks, 2L)

  exp_out <- readRDS("t1.l2.AdjPValues.rds")
  expected_cols_look2 <- colnames(exp_out)
  testthat::expect_true(all(expected_cols_look2 %in% colnames(state$mcpObj$AdjPValues)))
  testthat::expect_equal(state$mcpObj$AdjPValues, exp_out)

  exp_out <- readRDS("t1.l2.IndexSet.rds")
  testthat::expect_equal(state$mcpObj$IndexSet, exp_out)

  # Look 3
  state <- AnalyzeLook_PC(
    state,
    p_raw = c(H2 = 0.05, H4 = 0.10),
    plotGraphs = FALSE
  )

  testthat::expect_equal(state$completed_looks, 3L)

  exp_out <- readRDS("t1.l3.AdjPValues.rds")
  expected_cols_look3 <- colnames(exp_out)
  testthat::expect_true(all(expected_cols_look3 %in% colnames(state$mcpObj$AdjPValues)))
  testthat::expect_equal(state$mcpObj$AdjPValues, exp_out)


  # TODO: Debug this
  testthat::expect_true(state$trial_completed)

  # Plot helper should run
  g <- PlotAnalysisGraph(state, stage = 1)
  testthat::expect_true(!is.null(g))

  # print() should not error
  testthat::expect_no_error(print(state))
})

############
# Test2
testthat::test_that("Test 2: PC analysis API scaffolds (strategy modification)", {
  testthat::skip_if_not(
    isTRUE(getOption("AdaptGMCP.run_pc_analysis_api_tests", FALSE)),
    "Scaffold only: set options(AdaptGMCP.run_pc_analysis_api_tests = TRUE) and fill expected values"
  )
  wi <- c(1 / 2, 1 / 2, 0, 0)
  g <- matrix(c(
    0, 1 / 2, 1 / 2, 0,
    1 / 2, 0, 0, 1 / 2,
    0, 1, 0, 0,
    1, 0, 0, 0
  ), byrow = TRUE, nrow = 4)

  t <- c(0.5, 0.7, 1)
  alp <- 0.025

  corr <- matrix(c(
    1, 0.5, 0.5, NA,
    0.5, 1, NA, 0.5,
    0.5, NA, 1, 0.5,
    NA, 0.5, 0.5, 1
  ), byrow = TRUE, nrow = 4)

  tt <- "Partly-Parametric"
  des <- "asOF"

  #browser()

  state <- SetupAnalysis_PC(
    WI = wi,
    G = g,
    test.type = tt,
    alpha = alp,
    info_frac = t,
    typeOfDesign = des,
    Correlation = corr,
    plotGraphs = FALSE,
    MultipleWinners = FALSE
  )

  # Look 1
  state <- AnalyzeLook_PC(
    state,
    p_raw = c(H1 = 0.03, H2 = 0.20, H3 = 0.10, H4 = 0.25),
    plotGraphs = FALSE
  )

  exp_out <- readRDS("t2.l1.AdjPValues.rds")
  testthat::expect_equal(state$mcpObj$AdjPValues, exp_out)

  # Apply a strategy modification before analyzing look 2
  new_weights <- c(H1 = 0.5, H2 = 0.25, H3 = 0.25, H4 = 0)
  new_G <- matrix(c(
    0, 0.5, 0.5, 0,
    0.5, 0, 0.5, 0,
    0.5, 0.5, 0, 0,
    1, 0, 0, 0
  ), byrow = TRUE, nrow = 4)

  state <- AnalyzeLook_PC(
    state,
    p_raw = c(H1 = 0.02, H2 = 0.12, H3 = 0.08, H4 = 0.20),
    new_weights = new_weights,
    new_G = new_G,
    plotGraphs = FALSE
  )

  testthat::expect_equal(state$completed_looks, 2L)

  # TODO: Debug this
  exp_out <- readRDS("t2.l2.AdjPValues.rds")
  testthat::expect_equal(state$mcpObj$AdjPValues, exp_out)

  # TODO: Debug this
  # browser()
  exp_out <- readRDS("t2.l2.IndexSet.rds")
  testthat::expect_equal(state$mcpObj$IndexSet, exp_out)

  # Plot helper should run
  g <- PlotAnalysisGraph(state, stage = 2)
  testthat::expect_true(!is.null(g))
})

############
# Test3
testthat::test_that("Test 3: PC analysis API scaffolds (full transition at look 2)", {
  testthat::skip_if_not(
    isTRUE(getOption("AdaptGMCP.run_pc_analysis_api_tests", FALSE)),
    "Scaffold only: set options(AdaptGMCP.run_pc_analysis_api_tests = TRUE) and fill expected values"
  )

  wi <- c(1 / 2, 1 / 2, 0, 0)
  g <- matrix(c(
    0, 1 / 2, 1 / 2, 0,
    1 / 2, 0, 0, 1 / 2,
    0, 1, 0, 0,
    1, 0, 0, 0
  ), byrow = TRUE, nrow = 4)

  t <- c(0.5, 0.7, 1)
  alp <- 0.025

  corr <- matrix(c(
    1, 0.5, 0.5, NA,
    0.5, 1, NA, 0.5,
    0.5, NA, 1, 0.5,
    NA, 0.5, 0.5, 1
  ), byrow = TRUE, nrow = 4)

  tt <- "Partly-Parametric"
  des <- "asOF"

  state <- SetupAnalysis_PC(
    WI = wi,
    G = g,
    test.type = tt,
    alpha = alp,
    info_frac = t,
    typeOfDesign = des,
    Correlation = corr,
    plotGraphs = FALSE,
    MultipleWinners = TRUE
  )

  exp_out <- readRDS("t3.bdryTab.rds")
  testthat::expect_equal(state$mcpObj$bdryTab, exp_out)

  # Look 1
  state <- AnalyzeLook_PC(
    state,
    p_raw = c(H1 = 0.04, H2 = 0.18, H3 = 0.12, H4 = 0.22),
    plotGraphs = FALSE
  )

  exp_out <- readRDS("t3.l1.AdjPValues.rds")
  testthat::expect_equal(state$mcpObj$AdjPValues, exp_out)

  # Full transition at look 2:
  # 1) selection (drop H4)
  selection <- c("H1", "H2", "H3")

  # 2) strategy modification for continuing hypotheses
  new_weights <- c(H1 = 0.5, H2 = 0.5, H3 = 0)
  new_G <- matrix(c(
    0, 1, 0,
    0.5, 0, 0.5,
    1, 0, 0
  ), byrow = TRUE, nrow = 3)

  # 3) correlation update (example: tweak H1-H2 and H1-H3)
  new_correlation <- corr
  new_correlation[1, 2] <- new_correlation[2, 1] <- 0.3
  new_correlation[1, 3] <- new_correlation[3, 1] <- 0.4

  state <- AnalyzeLook_PC(
    state,
    p_raw = c(H1 = 0.03, H2 = 0.10, H3 = 0.08),
    selection = selection,
    new_weights = new_weights,
    new_G = new_G,
    new_correlation = new_correlation,
    plotGraphs = FALSE
  )

  testthat::expect_equal(state$completed_looks, 2L)

  exp_out <- readRDS("t3.l2.AdjPValues.rds")
  testthat::expect_equal(state$mcpObj$AdjPValues, exp_out)

  # TODO: Debug this
  # browser()
  g <- PlotAnalysisGraph(state, stage = 2)
  testthat::expect_true(!is.null(g))
})
