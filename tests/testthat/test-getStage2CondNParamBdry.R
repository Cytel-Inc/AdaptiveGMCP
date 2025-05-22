# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

test_that("Test Stage-2 adapted Non-Parametric boundary(Cum.)", {
  # Test Case-1: Consistency with East 2-sample Muller-Schefer Benchmarks
  a1 <- 0.00152532273778278
  p1 <- 1 - pnorm(1.6904003395167)
  v <- 1
  BJ <- 0.137062084206029 # CER/CP under Null
  SS1 <- 50
  SS2 <- 100
  benchmark <- 1 - pnorm(1.96859564104221)
  out <- getStage2CondNParamBdry(a1, p1, v, BJ, SS1, SS2)
  expect_equal(object = out$Stage2AdjBdry, expected = benchmark, tolerance = 1E-9)
  #----------------------------------------------------------------------------
  # Test Case-2: Consistency with East 2-sample Muller-Schefer Benchmarks
  a1 <- 0.00152532273778278
  p1 <- 1 - pnorm(1.75498559359717)
  v <- 1
  BJ <- 0.152051635710302 # CER/CP under Null
  SS1 <- 110
  SS2 <- 299
  benchmark <- 1 - pnorm(1.88152587781072)
  out <- getStage2CondNParamBdry(a1, p1, v, BJ, SS1, SS2)
  expect_equal(object = out$Stage2AdjBdry, expected = benchmark, tolerance = 1E-9)

  #----------------------------------------------------------------------------
  # Test Case-3: Consistency with planned boundary no adaptation is performed(R-pact benchmarks)
  alpha <- 0.025
  w <- c(0.2, 0.3, 0.4, 0.1)
  info_frac <- c(0.5, 1)
  PlanBdry <- lapply(1:length(w), function(x) {
    des <- rpact::getDesignGroupSequential(
      kMax = 2,
      alpha = alpha * w[x],
      informationRates = info_frac,
      typeOfDesign = "asOF"
    )
    des$stageLevels
  })

  Stage1PlanBdry <- unlist(lapply(1:length(PlanBdry), function(i) PlanBdry[[i]][1]))
  Stage2PlanBdry <- unlist(lapply(1:length(PlanBdry), function(i) PlanBdry[[i]][2]))

  ss1 <- c(133, 133, 133, 133)
  ss2 <- c(266, 266, 266, 266)
  p1 <- c(0.01, 0.02, 0.03, 0.04)
  pcer <- unlist(lapply(1:length(ss1), function(x) {
    getPCER(a2 = Stage2PlanBdry[x], p1 = p1[x], ss1 = ss1[x], ss2 = ss2[x])
  }))
  out <- getStage2CondNParamBdry(a1 = Stage1PlanBdry, p1 = p1, v = w, BJ = sum(pcer), SS1 = ss1, SS2 = ss2)
  expect_equal(object = out$Stage2AdjBdry, expected = Stage2PlanBdry)

  #------------------------------------------------------------------------------
  # Test Case-4: Consistency with planned boundary no adaptation is performed(R-pact benchmarks)
  alpha <- 0.025
  w <- c(0.2, 0.3, 0.4, 0.1)
  info_frac <- c(0.5, 1)
  PlanBdry <- lapply(1:length(w), function(x) {
    des <- rpact::getDesignGroupSequential(
      kMax = 2,
      alpha = alpha * w[x],
      informationRates = info_frac,
      typeOfDesign = "asOF"
    )
    des$stageLevels
  })

  Stage1PlanBdry <- unlist(lapply(1:length(PlanBdry), function(i) PlanBdry[[i]][1]))
  Stage2PlanBdry <- unlist(lapply(1:length(PlanBdry), function(i) PlanBdry[[i]][2]))

  ss1 <- c(133, 133, 133, 133)
  ss2 <- c(266, 266, 266, 266)
  p1 <- c(0.01, 0.02, 0.03, 0.04)
  pcer <- unlist(lapply(1:length(ss1), function(x) {
    getPCER(a2 = Stage2PlanBdry[x], p1 = p1[x], ss1 = ss1[x], ss2 = ss2[x])
  }))
  out <- getStage2CondNParamBdry(a1 = Stage1PlanBdry, p1 = p1, v = w, BJ = sum(pcer), SS1 = ss1, SS2 = ss2)
  expect_equal(object = out$Stage2AdjBdry, expected = Stage2PlanBdry)

  #------------------------------------------------------------------------------
  # Test Case-5: Consistency with the plan stage-2 boundary for No Early efficacy
  alpha <- 0.025
  w <- c(0.8, 0.2)
  info_frac <- c(0.5, 1)
  PlanBdry <- lapply(1:length(w), function(x) {
    des <- rpact::getDesignGroupSequential(
      kMax = 2,
      alpha = alpha * w[x],
      informationRates = info_frac,
      typeOfDesign = "noEarlyEfficacy"
    )
    des$stageLevels
  })

  a1 <- unlist(lapply(1:length(PlanBdry), function(i) PlanBdry[[i]][1]))
  a2 <- unlist(lapply(1:length(PlanBdry), function(i) PlanBdry[[i]][2]))

  p1 <- c(0.01, 0.02)
  SS1 <- c(50, 50)
  SS2 <- c(100, 100)
  pcer <- unlist(lapply(1:length(SS1), function(x) {
    getPCER(a2 = a2[x], p1 = p1[x], ss1 = SS1[x], ss2 = SS2[x])
  }))
  BJ <- sum(pcer)

  v <- c(0.8, 0.2)
  out <- getStage2CondNParamBdry(a1 = a1, p1 = p1, v = v, BJ = BJ, SS1 = SS1, SS2 = SS2)

  expect_equal(object = out$Stage2AdjBdry, expected = a2, tolerance = 1E-9)

  #------------------------------------------------------------------------------
  # Test Case-6: BJ >= 1 No Early efficacy with very small stage-1 p values)
  # (=> HJ could be rejected at level α on the basis of stage one data alone)

  alpha <- 0.025
  w <- c(0.8, 0.2)
  info_frac <- c(0.5, 1)
  PlanBdry <- lapply(1:length(w), function(x) {
    des <- rpact::getDesignGroupSequential(
      kMax = 2,
      alpha = alpha * w[x],
      informationRates = info_frac,
      typeOfDesign = "noEarlyEfficacy"
    )
    des$stageLevels
  })

  a1 <- unlist(lapply(1:length(PlanBdry), function(i) PlanBdry[[i]][1]))
  a2 <- unlist(lapply(1:length(PlanBdry), function(i) PlanBdry[[i]][2]))

  p1 <- c(0.0001, 0.002)
  SS1 <- c(50, 50)
  SS2 <- c(100, 100)
  pcer <- unlist(lapply(1:length(SS1), function(x) {
    getPCER(a2 = a2[x], p1 = p1[x], ss1 = SS1[x], ss2 = SS2[x])
  }))
  BJ <- sum(pcer)

  v <- c(0.8, 0.2)
  out <- getStage2CondNParamBdry(a1 = a1, p1 = p1, v = v, BJ = BJ, SS1 = SS1, SS2 = SS2)

  expect_equal(object = out$Stage2AdjBdry, expected = a2, tolerance = 1E-9)

  #-------------------------------------------------------------------------------
  # Test Case-7: Stage-2 boundaries are not greater than 1(Early Efficacy, BJ >= 1)
  alpha <- 0.025
  w <- c(0.9, 0.1)
  info_frac <- c(0.5, 1)
  PlanBdry <- lapply(1:length(w), function(x) {
    des <- rpact::getDesignGroupSequential(
      kMax = 2,
      alpha = alpha * w[x],
      informationRates = info_frac,
      typeOfDesign = "noEarlyEfficacy"
    )
    des$stageLevels
  })

  a1 <- unlist(lapply(1:length(PlanBdry), function(i) PlanBdry[[i]][1]))
  a2 <- unlist(lapply(1:length(PlanBdry), function(i) PlanBdry[[i]][2]))

  p1 <- c(0.0000001, 0.00000002)
  SS1 <- c(50, 50)
  SS2 <- c(100, 100)
  pcer <- unlist(lapply(1:length(SS1), function(x) {
    getPCER(a2 = a2[x], p1 = p1[x], ss1 = SS1[x], ss2 = SS2[x])
  }))
  BJ <- sum(pcer)

  v <- c(0.1, 0.9)
  out <- getStage2CondNParamBdry(a1 = a1, p1 = p1, v = v, BJ = BJ, SS1 = SS1, SS2 = SS2)
  expect_equal(object = out$Stage2AdjBdry <= 1, expected = c(T, T))
})
