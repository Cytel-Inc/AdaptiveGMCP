# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

test_that("Test genNornToOther", {
  #----------------------------------------------------
  # Test Case-4: Mixed(both Continuous)
  nArmID <- 3
  nSubject <- 100000
  vEPs <- c(1, 3)
  vEPType <- c("Continuous", "Continuous")
  lNormMean <- list(
    "EP1" = c(0, 0.4, 0.3),
    "EP2" = NA,
    "Ep3" = c(0.2, 0.55, 0.35),
    "EP4" = NA
  )
  lNormStdDev <- list(
    "EP1" = c(1.1, 1.2, 1.3),
    "EP2" = NA,
    "Ep3" = c(1.15, 1.25, 1.35),
    "EP4" = NA
  )
  lProp <- list(
    "EP1" = NA,
    "EP2" = c(0.1, 0.2, 0.3),
    "EP3" = NA,
    "EP4" = c(0.5, 0.6, 0.7)
  )
  mNormCorr <- matrix(
    c(
      1, 0.5, 0.5, 0.5,
      0.5, 1, 0.5, 0.5,
      0.5, 0.5, 1, 0.5,
      0.5, 0.5, 0.5, 1
    ),
    nrow = 4
  )
  nSeed <- 123

  out <- genNormToOther2(
    nArmID = nArmID, nSubject = nSubject, vEPs = vEPs,
    vEPType = vEPType, lNormMean = lNormMean, lNormStdDev = lNormStdDev,
    lProp = lProp, mNormCorr = mNormCorr, nSeed = nSeed
  )

  expect_equal(
    object = colSums(out) / nSubject - c(0.3, 0.35),
    expected = c(0, 0),
    tolerance = 1E-2
  )
  #------------------------------------------------------------------
  # Test Case-5: Mixed(Single Continuous)
  vEPs <- c(1)
  vEPType <- c("Continuous")
  out <- genNormToOther2(
    nArmID = nArmID, nSubject = nSubject, vEPs = vEPs,
    vEPType = vEPType, lNormMean = lNormMean, lNormStdDev = lNormStdDev,
    lProp = lProp, mNormCorr = mNormCorr, nSeed = nSeed
  )

  expect_equal(
    object = colSums(out) / nSubject - c(0.3),
    expected = c(0),
    tolerance = 1E-2
  )
  #-----------------------------------------------------------------
  # Test Case-6: Mixed(Single Binary)
  vEPs <- c(2)
  vEPType <- c("Binary")
  out <- genNormToOther2(
    nArmID = nArmID, nSubject = nSubject, vEPs = vEPs,
    vEPType = vEPType, lNormMean = lNormMean, lNormStdDev = lNormStdDev,
    lProp = lProp, mNormCorr = mNormCorr, nSeed = nSeed
  )

  expect_equal(
    object = colSums(out) / nSubject - c(0.3),
    expected = c(0),
    tolerance = 1E-2
  )
  #-----------------------------------------------------------------
  # Test Case-7: Mixed(Both Binary)
  vEPs <- c(2, 4)
  vEPType <- c("Binary", "Binary")
  out <- genNormToOther2(
    nArmID = nArmID, nSubject = nSubject, vEPs = vEPs,
    vEPType = vEPType, lNormMean = lNormMean, lNormStdDev = lNormStdDev,
    lProp = lProp, mNormCorr = mNormCorr, nSeed = nSeed
  )
  expect_equal(
    object = colSums(out) / nSubject - c(0.3, 0.7),
    expected = c(0, 0),
    tolerance = 1E-2
  )
  #----------------------------------------------------------------
  # Test Case-4: Mixed(Continuous Binary)
  vEPs <- c(1, 3, 2, 4)
  vEPType <- c("Continuous", "Continuous", "Binary", "Binary")
  out <- genNormToOther2(
    nArmID = nArmID, nSubject = nSubject, vEPs = vEPs,
    vEPType = vEPType, lNormMean = lNormMean, lNormStdDev = lNormStdDev,
    lProp = lProp, mNormCorr = mNormCorr, nSeed = nSeed
  )

  expect_equal(
    object = colSums(out) / nSubject - c(0.3, 0.35, 0.3, 0.7),
    expected = c(0, 0, 0, 0),
    tolerance = 1E-2
  )
  #---------------------------------------------------------------
})
