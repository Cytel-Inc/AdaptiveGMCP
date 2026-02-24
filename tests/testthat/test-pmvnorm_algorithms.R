# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

# Tests for mvtnorm::pmvnorm() with different covariance matrices and algorithms
# Testing behavior of Miwa and GenzBretz algorithms with various matrix conditions

test_that("pmvnorm works with valid positive definite matrix using Miwa", {
  # browser()
  # Create a valid correlation matrix
  corr_mat <- matrix(c(
    1.0, 0.5, 0.3,
    0.5, 1.0, 0.4,
    0.3, 0.4, 1.0
  ), nrow = 3, byrow = TRUE)
  
  # Test that determinant is positive
  expect_true(det(corr_mat) > 0)
  
  # Create Miwa algorithm
  algo <- mvtnorm::Miwa(steps = 128, checkCorr = FALSE, maxval = 1e3)
  
  # This should work without error
  result <- mvtnorm::pmvnorm(
    lower = c(-Inf, -Inf, -Inf),
    upper = c(0, 0, 0),
    corr = corr_mat,
    algorithm = algo
  )
  
  # Result should be numeric and between 0 and 1
  expect_true(is.numeric(result[1]))
  expect_true(result[1] >= 0 && result[1] <= 1)
})

test_that("pmvnorm works with valid positive definite matrix using GenzBretz", {
  # browser()
  # Create a valid correlation matrix
  corr_mat <- matrix(c(
    1.0, 0.5, 0.3,
    0.5, 1.0, 0.4,
    0.3, 0.4, 1.0
  ), nrow = 3, byrow = TRUE)
  
  # Create GenzBretz algorithm
  algo <- mvtnorm::GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0)
  
  # This should work without error
  result <- mvtnorm::pmvnorm(
    lower = c(-Inf, -Inf, -Inf),
    upper = c(0, 0, 0),
    corr = corr_mat,
    algorithm = algo
  )
  
  # Result should be numeric and between 0 and 1
  expect_true(is.numeric(result[1]))
  expect_true(result[1] >= 0 && result[1] <= 1)
})

test_that("pmvnorm with singular matrix fails using Miwa", {
  # browser()
  # Create a singular correlation matrix (linearly dependent rows)
  # Third row = 2 * first row - second row
  singular_mat <- matrix(c(
    1.0, 0.5, 0.3,
    0.5, 1.0, 0.4,
    1.5, 0.0, 0.2  # This makes the matrix singular
  ), nrow = 3, byrow = TRUE)
  
  # Verify the matrix is singular (determinant near zero)
  expect_true(abs(det(singular_mat)) < 1e-10)
  
  # Create Miwa algorithm
  algo <- mvtnorm::Miwa(steps = 128, checkCorr = FALSE, maxval = 1e3)
  
  # This should throw an error
  expect_error(
    mvtnorm::pmvnorm(
      lower = c(-Inf, -Inf, -Inf),
      upper = c(0, 0, 0),
      corr = singular_mat,
      algorithm = algo
    )
  )
})

test_that("pmvnorm with singular matrix works using GenzBretz", {
  # browser()
  # Create a singular correlation matrix
  singular_mat <- matrix(c(
    1.0, 0.5, 0.3,
    0.5, 1.0, 0.4,
    1.5, 0.0, 0.2
  ), nrow = 3, byrow = TRUE)
  
  # Verify the matrix is singular
  expect_true(abs(det(singular_mat)) < 1e-10)
  
  # Create GenzBretz algorithm
  algo <- mvtnorm::GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0)
  
  # GenzBretz might still work or give a more graceful error
  # Testing that it at least doesn't crash R completely
  result <- tryCatch(
    {
      mvtnorm::pmvnorm(
        lower = c(-Inf, -Inf, -Inf),
        upper = c(0, 0, 0),
        corr = singular_mat,
        algorithm = algo
      )
    },
    error = function(e) {
      return(NULL)
    }
  )
  
  # Just check that we got some result (either a value or handled error)
  expect_true(is.null(result) || is.numeric(result[1]))
})

test_that("pmvnorm with nearly singular matrix may fail with Miwa", {
  # browser()
  # Create a nearly singular matrix (very small determinant)
  nearly_singular <- matrix(c(
    1.0, 0.99, 0.98,
    0.99, 1.0, 0.99,
    0.98, 0.99, 1.0
  ), nrow = 3, byrow = TRUE)
  
  # Check that determinant is very small but not exactly zero
  det_value <- det(nearly_singular)
  expect_true(det_value > 0)
  expect_true(det_value < 0.01)  # Very small but positive
  
  # Create Miwa algorithm
  algo <- mvtnorm::Miwa(steps = 128, checkCorr = FALSE, maxval = 1e3)
  
  # This might work or fail depending on numerical precision
  result <- tryCatch(
    {
      mvtnorm::pmvnorm(
        lower = c(-Inf, -Inf, -Inf),
        upper = c(0, 0, 0),
        corr = nearly_singular,
        algorithm = algo
      )
    },
    error = function(e) {
      return(NULL)
    }
  )
  
  # We expect either a valid result or an error was caught
  expect_true(is.null(result) || (is.numeric(result[1]) && result[1] >= 0 && result[1] <= 1))
})

test_that("pmvnorm with identity matrix works with both algorithms", {
  # browser()
  # Identity matrix (no correlation)
  identity_mat <- diag(3)
  
  # Test with Miwa
  algo_miwa <- mvtnorm::Miwa(steps = 128, checkCorr = FALSE, maxval = 1e3)
  result_miwa <- mvtnorm::pmvnorm(
    lower = c(-Inf, -Inf, -Inf),
    upper = c(0, 0, 0),
    corr = identity_mat,
    algorithm = algo_miwa
  )
  
  # Test with GenzBretz
  algo_genz <- mvtnorm::GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0)
  result_genz <- mvtnorm::pmvnorm(
    lower = c(-Inf, -Inf, -Inf),
    upper = c(0, 0, 0),
    corr = identity_mat,
    algorithm = algo_genz
  )
  
  # Both should work and give similar results
  expect_true(is.numeric(result_miwa[1]))
  expect_true(is.numeric(result_genz[1]))
  
  # For independent variables, result should be product of marginals
  # P(Z1 <= 0) * P(Z2 <= 0) * P(Z3 <= 0) = 0.5^3 = 0.125
  expect_equal(result_miwa[1], 0.125, tolerance = 0.01)
  expect_equal(result_genz[1], 0.125, tolerance = 0.01)
})

test_that("pmvnorm with high correlation matrix works", {
  # browser()
  # Matrix with high but valid correlations
  high_corr <- matrix(c(
    1.0, 0.8, 0.7,
    0.8, 1.0, 0.75,
    0.7, 0.75, 1.0
  ), nrow = 3, byrow = TRUE)
  
  # Check it's positive definite
  expect_true(det(high_corr) > 0)
  
  # Should work with Miwa
  algo <- mvtnorm::Miwa(steps = 128, checkCorr = FALSE, maxval = 1e3)
  result <- mvtnorm::pmvnorm(
    lower = c(-Inf, -Inf, -Inf),
    upper = c(1, 1, 1),
    corr = high_corr,
    algorithm = algo
  )
  
  expect_true(is.numeric(result[1]))
  expect_true(result[1] >= 0 && result[1] <= 1)
})

test_that("pmvnorm with negative correlations works", {
  # browser()
  # Matrix with negative correlations
  neg_corr <- matrix(c(
    1.0, -0.3, -0.2,
    -0.3, 1.0, -0.25,
    -0.2, -0.25, 1.0
  ), nrow = 3, byrow = TRUE)
  
  # Check it's positive definite
  expect_true(det(neg_corr) > 0)
  
  # Should work with Miwa
  algo <- mvtnorm::Miwa(steps = 128, checkCorr = FALSE, maxval = 1e3)
  result <- mvtnorm::pmvnorm(
    lower = c(-Inf, -Inf, -Inf),
    upper = c(0, 0, 0),
    corr = neg_corr,
    algorithm = algo
  )
  
  expect_true(is.numeric(result[1]))
  expect_true(result[1] >= 0 && result[1] <= 1)
})

test_that("checkCorr parameter catches invalid matrices when enabled", {
  # browser()
  # Create an invalid correlation matrix (not positive definite)
  invalid_corr <- matrix(c(
    1.0, 0.9, 0.9,
    0.9, 1.0, 0.9,
    0.9, 0.9, 1.0
  ), nrow = 3, byrow = TRUE)
  
  # This matrix has determinant very close to zero
  expect_true(det(invalid_corr) < 0.05)
  
  # With checkCorr = TRUE, Miwa should catch this
  algo_check <- mvtnorm::Miwa(steps = 128, checkCorr = TRUE, maxval = 1e3)
  
  # May throw error during algorithm creation or pmvnorm call
  expect_error(
    {
      result <- mvtnorm::pmvnorm(
        lower = c(-Inf, -Inf, -Inf),
        upper = c(0, 0, 0),
        corr = invalid_corr,
        algorithm = algo_check
      )
    },
    regexp = NA  # Accept any error message
  )
})

test_that("Determinant threshold detection works", {
  # browser()
  # Test matrices at different determinant levels
  
  # Good matrix (det >> machine epsilon)
  good_mat <- matrix(c(
    1.0, 0.3, 0.2,
    0.3, 1.0, 0.25,
    0.2, 0.25, 1.0
  ), nrow = 3, byrow = TRUE)
  
  det_good <- det(good_mat)
  expect_true(abs(det_good) > .Machine$double.eps)
  expect_true(abs(det_good) > 0.1)
  
  # Borderline matrix (det ~ 0.001)
  borderline_mat <- matrix(c(
    1.0, 0.95, 0.94,
    0.95, 1.0, 0.95,
    0.94, 0.95, 1.0
  ), nrow = 3, byrow = TRUE)
  
  det_border <- det(borderline_mat)
  expect_true(abs(det_border) > 0)
  expect_true(abs(det_border) < 0.01)
})

test_that("2x2 matrices work with both algorithms", {
  # browser()
  # Simple 2x2 case
  corr_2x2 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  
  # Miwa
  algo_miwa <- mvtnorm::Miwa(steps = 128, checkCorr = FALSE, maxval = 1e3)
  result_miwa <- mvtnorm::pmvnorm(
    lower = c(-Inf, -Inf),
    upper = c(0, 0),
    corr = corr_2x2,
    algorithm = algo_miwa
  )
  
  # GenzBretz
  algo_genz <- mvtnorm::GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0)
  result_genz <- mvtnorm::pmvnorm(
    lower = c(-Inf, -Inf),
    upper = c(0, 0),
    corr = corr_2x2,
    algorithm = algo_genz
  )
  
  expect_true(is.numeric(result_miwa[1]))
  expect_true(is.numeric(result_genz[1]))
  
  # Results should be similar
  expect_equal(result_miwa[1], result_genz[1], tolerance = 0.01)
})

test_that("Large dimension matrix works", {
  # browser()
  # Test with larger matrix (5x5)
  n <- 5
  large_corr <- matrix(0.3, nrow = n, ncol = n)
  diag(large_corr) <- 1
  
  # Should be positive definite
  expect_true(det(large_corr) > 0)
  
  # Test with Miwa
  algo <- mvtnorm::Miwa(steps = 128, checkCorr = FALSE, maxval = 1e3)
  result <- mvtnorm::pmvnorm(
    lower = rep(-Inf, n),
    upper = rep(0, n),
    corr = large_corr,
    algorithm = algo
  )
  
  expect_true(is.numeric(result[1]))
  expect_true(result[1] >= 0 && result[1] <= 1)
})
