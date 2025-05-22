# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

test_that("Test Non-Parametric Boundary Computation", {
  # Test Case-1: Spending Function : LD-OF
  nLooks <- 2
  sig_level <- 0.025
  info_frac <- c(1 / 2, 1)
  typeOfDesign <- "asOF"
  benchmark <- c(0.00152532273778278, 0.0244997717772641)
  outBdry <- unlist(getPlanNonParmBdry(nLooks = nLooks, sig_level = sig_level, info_frac = info_frac, typeOfDesign = typeOfDesign))
  names(outBdry) <- NULL
  expect_equal(object = outBdry, expected = benchmark, tolerance = 1E-6)

  # Test Case-2: sig_level = 0
  nLooks <- 2
  sig_level <- 0
  info_frac <- c(1 / 2, 1)
  typeOfDesign <- "asOF"
  benchmark <- c(0, 0)
  outBdry <- unlist(getPlanNonParmBdry(nLooks = nLooks, sig_level = sig_level, info_frac = info_frac, typeOfDesign = typeOfDesign))
  names(outBdry) <- NULL
  expect_equal(object = outBdry, expected = benchmark)

  # Test Case-3: no early efficacy
  nLooks <- 2
  sig_level <- 0.01
  info_frac <- c(1 / 2, 1)
  typeOfDesign <- "noEarlyEfficacy"
  benchmark <- c(0, sig_level)
  outBdry <- unlist(getPlanNonParmBdry(nLooks = nLooks, sig_level = sig_level, info_frac = info_frac, typeOfDesign = typeOfDesign))
  names(outBdry) <- NULL
  expect_equal(object = outBdry, expected = benchmark)

  # Test Case-4: Negative Case,with No Info Frac Error is expected
  nLooks <- 2
  sig_level <- 0.01
  info_frac <- c(0, 1)
  typeOfDesign <- "asOF"
  expect_error(getPlanNonParmBdry(nLooks = nLooks, sig_level = sig_level, info_frac = info_frac, typeOfDesign = typeOfDesign))

  # Test Case-5: FSD
  nLooks <- 1
  sig_level <- 0.0001
  info_frac <- 1
  typeOfDesign <- "asOF"
  outBdry <- unlist(getPlanNonParmBdry(nLooks = nLooks, sig_level = sig_level, info_frac = info_frac, typeOfDesign = typeOfDesign))
  names(outBdry) <- NULL
  expect_equal(object = outBdry[1], expected = sig_level)
})
