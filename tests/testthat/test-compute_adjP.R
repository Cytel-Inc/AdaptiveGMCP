# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

test_that("Test computations of Adjusted P-values", {
  # Test Case-1, Parametric Dunnett's test

  test.type <- "Dunnett"
  p <- matrix(
    c(
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.649907459, 0.087058602, 0.261696810, 0.00053251,
      0.649907459, 0.087058602, 0.261696810, 0.00053251,
      0.649907459, 0.087058602, 0.261696810, 0.00053251,
      0.649907459, 0.087058602, 0.261696810, 0.00053251
    ),
    nrow = 8, byrow = T
  )

  h <- matrix(
    c(
      1, 1, 1, 1, 0.5, 0.5, 0, 0,
      1, 1, 1, 0, 0.75, 0, 0.25, 0,
      1, 0, 1, 0, 0.5, 0, 0.5, 0,
      0, 0, 0, 1, 0, 0, 0, 1,
      1, 1, 1, 1, 1 / 4, 1 / 4, 1 / 4, 1 / 4,
      1, 1, 1, 0, 1 / 3, 1 / 3, 1 / 3, 0,
      1, 0, 1, 0, 1 / 2, 0, 1 / 2, 0,
      0, 0, 0, 1, 0, 0, 0, 1
    ),
    nrow = 8, byrow = T
  )
  cr <- matrix(0.5, nrow = 4, ncol = 4)
  diag(cr) <- 1

  adjP <- c()
  for (i in 1:nrow(p))
  {
    adjP[i] <- compute_adjP(h = h[i, ], cr = cr, p = p[i, ], test.type = test.type)$adj_pj
  }

  OutRej <- adjP < c(rep(0.002, 4), rep(0.025, 4))
  benchmark <- c(F, F, F, F, T, F, F, T)
  expect_equal(object = OutRej, expected = benchmark)

  #-----------------------------------------------------------------
  # Test Case-2, Partly-Parametric test
  test.type <- "Partly-Parametric"
  p <- matrix(
    c(
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.649907459, 0.087058602, 0.261696810, 0.00053251,
      0.649907459, 0.087058602, 0.261696810, 0.00053251,
      0.649907459, 0.087058602, 0.261696810, 0.00053251,
      0.649907459, 0.087058602, 0.261696810, 0.00053251
    ),
    nrow = 8, byrow = T
  )

  h <- matrix(
    c(
      1, 1, 1, 1, 0.5, 0.5, 0, 0,
      1, 1, 1, 0, 0.75, 0, 0.25, 0,
      1, 0, 1, 0, 0.5, 0, 0.5, 0,
      0, 0, 0, 1, 0, 0, 0, 1,
      1, 1, 1, 1, 1 / 4, 1 / 4, 1 / 4, 1 / 4,
      1, 1, 1, 0, 1 / 3, 1 / 3, 1 / 3, 0,
      1, 0, 1, 0, 1 / 2, 0, 1 / 2, 0,
      0, 0, 0, 1, 0, 0, 0, 1
    ),
    nrow = 8, byrow = T
  )
  cr <- matrix(c(
    1, 0.5, NA, NA,
    0.5, 1, NA, NA,
    NA, NA, 1, 0.5,
    NA, NA, 0.5, 1
  ), nrow = 4, byrow = T)

  adjP <- c()
  for (i in 1:nrow(p))
  {
    adjP[i] <- compute_adjP(h = h[i, ], cr = cr, p = p[i, ], test.type = test.type)$adj_pj
  }

  OutRej <- adjP < c(rep(0.002, 4), rep(0.025, 4))
  benchmark <- c(F, F, F, F, T, F, F, T)
  expect_equal(object = OutRej, expected = benchmark)

  #-----------------------------------------------------------------
  # Test Case-3, Bonferroni test
  test.type <- "Bonf"
  p <- matrix(
    c(
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.649907459, 0.087058602, 0.261696810, 0.00053251,
      0.649907459, 0.087058602, 0.261696810, 0.00053251,
      0.649907459, 0.087058602, 0.261696810, 0.00053251,
      0.649907459, 0.087058602, 0.261696810, 0.00053251
    ),
    nrow = 8, byrow = T
  )

  h <- matrix(
    c(
      1, 1, 1, 1, 0.5, 0.5, 0, 0,
      1, 1, 1, 0, 0.75, 0, 0.25, 0,
      1, 0, 1, 0, 0.5, 0, 0.5, 0,
      0, 0, 0, 1, 0, 0, 0, 1,
      1, 1, 1, 1, 1 / 4, 1 / 4, 1 / 4, 1 / 4,
      1, 1, 1, 0, 1 / 3, 1 / 3, 1 / 3, 0,
      1, 0, 1, 0, 1 / 2, 0, 1 / 2, 0,
      0, 0, 0, 1, 0, 0, 0, 1
    ),
    nrow = 8, byrow = T
  )
  cr <- matrix(c(
    1, NA, NA, NA,
    NA, 1, NA, NA,
    NA, NA, 1, NA,
    NA, NA, NA, 1
  ), nrow = 4, byrow = T)

  adjP <- c()
  for (i in 1:nrow(p))
  {
    adjP[i] <- compute_adjP(h = h[i, ], cr = cr, p = p[i, ], test.type = test.type)$adj_pj
  }

  OutRej <- adjP < c(rep(0.002, 4), rep(0.025, 4))
  benchmark <- c(F, F, F, F, T, F, F, T)
  expect_equal(object = OutRej, expected = benchmark)

  #-----------------------------------------------------------------
  # Test Case-4, Sidak test
  test.type <- "Sidak"
  p <- matrix(
    c(
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.649907459, 0.087058602, 0.261696810, 0.00053251,
      0.649907459, 0.087058602, 0.261696810, 0.00053251,
      0.649907459, 0.087058602, 0.261696810, 0.00053251,
      0.649907459, 0.087058602, 0.261696810, 0.00053251
    ),
    nrow = 8, byrow = T
  )

  h <- matrix(
    c(
      1, 1, 1, 1, 0.5, 0.5, 0, 0,
      1, 1, 1, 0, 0.75, 0, 0.25, 0,
      1, 0, 1, 0, 0.5, 0, 0.5, 0,
      0, 0, 0, 1, 0, 0, 0, 1,
      1, 1, 1, 1, 1 / 4, 1 / 4, 1 / 4, 1 / 4,
      1, 1, 1, 0, 1 / 3, 1 / 3, 1 / 3, 0,
      1, 0, 1, 0, 1 / 2, 0, 1 / 2, 0,
      0, 0, 0, 1, 0, 0, 0, 1
    ),
    nrow = 8, byrow = T
  )

  adjP <- c()
  for (i in 1:nrow(p))
  {
    adjP[i] <- compute_adjP(h = h[i, ], p = p[i, ], test.type = test.type)$adj_pj
  }

  OutRej <- adjP < c(rep(0.002, 4), rep(0.025, 4))
  benchmark <- c(F, F, F, F, T, F, F, T)
  expect_equal(object = OutRej, expected = benchmark)

  # Test Case-4, Sidak test
  test.type <- "Simes"
  p <- matrix(
    c(
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.252463317, 0.238189947, 0.071954653, 0.00361405,
      0.649907459, 0.087058602, 0.261696810, 0.00053251,
      0.649907459, 0.087058602, 0.261696810, 0.00053251,
      0.649907459, 0.087058602, 0.261696810, 0.00053251,
      0.649907459, 0.087058602, 0.261696810, 0.00053251
    ),
    nrow = 8, byrow = T
  )

  h <- matrix(
    c(
      1, 1, 1, 1, 0.5, 0.5, 0, 0,
      1, 1, 1, 0, 0.75, 0, 0.25, 0,
      1, 0, 1, 0, 0.5, 0, 0.5, 0,
      0, 0, 0, 1, 0, 0, 0, 1,
      1, 1, 1, 1, 1 / 4, 1 / 4, 1 / 4, 1 / 4,
      1, 1, 1, 0, 1 / 3, 1 / 3, 1 / 3, 0,
      1, 0, 1, 0, 1 / 2, 0, 1 / 2, 0,
      0, 0, 0, 1, 0, 0, 0, 1
    ),
    nrow = 8, byrow = T
  )

  adjP <- c()
  for (i in 1:nrow(p))
  {
    adjP[i] <- compute_adjP(h = h[i, ], p = p[i, ], test.type = test.type)$adj_pj
  }

  OutRej <- adjP < c(rep(0.002, 4), rep(0.025, 4))
  benchmark <- c(F, F, F, F, T, F, F, T)
  expect_equal(object = OutRej, expected = benchmark)
  #------------------------------------------------------------------------------
})
