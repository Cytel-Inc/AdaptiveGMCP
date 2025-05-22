# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

test_that("Test Inverse Normal Weights", {
  planSSIncr <- matrix(c(
    12, 15, 19,
    10, 20, 15
  ), nrow = 2, byrow = T)
  invNorm <- getInvNormWeights(planSSIncr = planSSIncr)
  expect_equal(object = sum(invNorm$W_Norm^2), expected = 1, tolerance = 1E-6)

  planSSIncr <- matrix(c(
    12, 15, 19,
    10, 20, 15,
    15, 15, 20
  ), nrow = 3, byrow = T)
  invNorm <- getInvNormWeights(planSSIncr = planSSIncr)
  expect_equal(object = rowSums(invNorm$W_Norm^2, na.rm = T), expected = c(1, 1), tolerance = 1E-6)
})




test_that("Test Combining p-value computations", {
  # Test Case-1: Independent Check for two Stage design
  planSSIncr <- matrix(c(
    10, 15, 19,
    80, 20, 15
  ), nrow = 2, byrow = T)
  invNorm <- getInvNormWeights(planSSIncr = planSSIncr)
  W_Norm <- invNorm$W_Norm
  adjPValue <- data.frame(matrix(c(0.01, 0.02), nrow = 1))
  colnames(adjPValue) <- rep("PAdj", 2)

  out <- CombinedPvalue(CurrentLook = 2, adjPValue = adjPValue, W_Norm = W_Norm)
  benchmark <- 1 - pnorm(sum(W_Norm * qnorm(1 - unlist(adjPValue))))

  expect_equal(object = out, expected = benchmark)

  # Test Case-2: Expect 0
  planSSIncr <- matrix(c(
    100000, 150000, 190000,
    150000, 190000, 10000
  ), nrow = 2, byrow = T)
  invNorm <- getInvNormWeights(planSSIncr = planSSIncr)
  W_Norm <- invNorm$W_Norm
  adjPValue <- data.frame(matrix(c(0.0, 0.0), nrow = 1))
  colnames(adjPValue) <- rep("PAdj", 2)
  out <- CombinedPvalue(CurrentLook = 2, adjPValue = adjPValue, W_Norm = W_Norm)
  expect_equal(object = out, expected = 0)

  # Test Case-2: Expect NA if one of the adjusted p-value is NA
  planSSIncr <- matrix(c(
    100000, 150000, 190000,
    150000, 190000, 10000
  ), nrow = 2, byrow = T)
  invNorm <- getInvNormWeights(planSSIncr = planSSIncr)
  W_Norm <- invNorm$W_Norm
  adjPValue <- data.frame(matrix(c(NA, 0.1), nrow = 1))
  colnames(adjPValue) <- rep("PAdj", 2)
  expect_equal(is.na(CombinedPvalue(CurrentLook = 2, adjPValue = adjPValue, W_Norm = W_Norm)),
    expected = TRUE
  )
})
