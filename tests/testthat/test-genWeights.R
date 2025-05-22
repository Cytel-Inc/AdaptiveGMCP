# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

test_that("Test Intersection Weights Computation", {
  ## Test Case-1 gate keeping strategy
  w <- c(1 / 2, 1 / 2, 0, 0)
  g <- matrix(c(
    0, 1 / 2, 1 / 2, 0,
    1 / 2, 0, 0, 1 / 2,
    0, 1, 0, 0,
    1, 0, 1, 0
  ), nrow = 4, byrow = T)

  benchmark <- gMCPLite::generateWeights(g = g, w = w)
  benchmark_idx <- as.vector(apply(benchmark[, 1:4], 1, function(x) {
    paste(x, collapse = "")
  }))

  out <- genWeights(w = w, g = g)
  outWeights <- out$IntersectionWeights
  rownames(outWeights) <- colnames(outWeights) <- NULL
  outWeights_idx <- as.vector(apply(outWeights[, 1:4], 1, function(x) {
    paste(x, collapse = "")
  }))
  new_idx <- unlist(lapply(outWeights_idx, function(x) which(x == benchmark_idx)))
  benchmark2 <- as.data.frame(benchmark[new_idx, ])
  rownames(benchmark2) <- colnames(benchmark2) <- NULL

  expect_equal(object = outWeights, expected = benchmark2)
  #----------------------------------------------------------
  ## Test Case-2 sequential strategy
  w <- c(1 / 4, 1 / 4, 1 / 4, 1 / 4)
  g <- matrix(c(
    0, 1 / 3, 1 / 3, 1 / 3,
    1 / 3, 0, 1 / 3, 1 / 3,
    1 / 3, 1 / 3, 0, 1 / 3,
    1 / 3, 1 / 3, 1 / 3, 0
  ), nrow = 4, byrow = T)

  benchmark <- gMCPLite::generateWeights(g = g, w = w)
  benchmark_idx <- as.vector(apply(benchmark[, 1:4], 1, function(x) {
    paste(x, collapse = "")
  }))

  out <- genWeights(w = w, g = g)
  outWeights <- out$IntersectionWeights
  rownames(outWeights) <- colnames(outWeights) <- NULL
  outWeights_idx <- as.vector(apply(outWeights[, 1:4], 1, function(x) {
    paste(x, collapse = "")
  }))
  new_idx <- unlist(lapply(outWeights_idx, function(x) which(x == benchmark_idx)))
  benchmark2 <- as.data.frame(benchmark[new_idx, ])
  rownames(benchmark2) <- colnames(benchmark2) <- NULL

  expect_equal(object = outWeights, expected = benchmark2)

  #----------------------------------------------------------
  ## Test Case-3 Error if dimension of w and g are not consistent
  w <- c(1 / 3, 1 / 3, 1 / 3)
  g <- matrix(c(
    0, 1,
    1, 0
  ), nrow = 2)
  expect_error(genWeights(w = w, g = g))

  w <- c(1 / 2, 1 / 2)
  g <- matrix(c(
    0, 1 / 2, 1 / 2,
    1 / 2, 0, 1 / 2,
    1 / 2, 1 / 2, 0
  ), nrow = 3)
  expect_error(genWeights(w = w, g = g))

  w <- c(1 / 3, 1 / 3, 1 / 3)
  g <- matrix(c(
    0, 1 / 2, 1 / 2,
    1 / 2, 0, 1 / 2,
    1 / 2, 1 / 2, 0
  ), nrow = 3)
  expect_error(genWeights(w = w, g = g, HypothesisName = c("A", "B")))
})
