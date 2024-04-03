test_that("Test Correlation Computation for Combining p-values(Dunnett)", {
  # Test Case-1: Consistency with getSigma results
  # Arms: 3, Eps: 1, Hypothesis: 2
  # getSigma Results
  sigma <- list("EP1" = c(1.5, 1.4, 1.6))
  allocRatio <- c(1, 2, 3)
  ctrSS <- c(133, 347)
  SS_Cum <- cbind(ctrSS, allocRatio[2] * ctrSS, allocRatio[3] * ctrSS)
  SS_Incr <- rbind(SS_Cum[1, ], SS_Cum[2, ] - SS_Cum[1, ])

  outgetSigma <- getSigma(SS_Cum, sigma, allocRatio)
  SigmaZ_Stage <- matrix(outgetSigma$SigmaZ$EP1, nrow = 4)
  sigmaZ_Stage1 <- SigmaZ_Stage[c(1, 2), c(1, 2)]
  sigmaZ_Stage2 <- SigmaZ_Stage[c(3, 4), c(3, 4)]

  # getPlanCorrelation Results
  outgetPlanCorr <- getPlanCorrelation(
    nHypothesis = 2,
    SS_Incr = SS_Incr,
    Arms.std.dev = sigma,
    test.type = "Dunnett"
  )

  corr_stage1 <- matrix(outgetPlanCorr$Stage1, nrow = 2)
  corr_stage2 <- matrix(outgetPlanCorr$Stage2, nrow = 2)

  expect_equal(object = corr_stage1, expected = sigmaZ_Stage1)
  expect_equal(object = corr_stage2, expected = sigmaZ_Stage2)
})
