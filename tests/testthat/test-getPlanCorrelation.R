test_that("Test Correlation Computation for Combining p-values(Dunnett)", {
  # Test Case-1: Consistency with getSigma results
  # Arms: 3, Eps: 1, Hypothesis: 2
  # getSigma Results
  CommonStdDev <- FALSE
  EpType <- list('EP1'='Continuous')
  sigma <- list("EP1" = c(1.5, 1.4, 1.6))
  prop.ctr <- NA
  allocRatio <- c(1, 2, 3)
  ctrSS <- c(133, 347)
  SS_Cum <- data.frame('ctr' = ctrSS,
                       'trt1' = allocRatio[2] * ctrSS,
                       'trt2' = allocRatio[3] * ctrSS)
  SS_Incr <- rbind(SS_Cum[1, ], SS_Cum[2, ] - SS_Cum[1, ])

  outgetSigma <- getSigma(EpType = EpType,
                          SS_Cum = SS_Cum,
                          sigma = sigma,
                          prop.ctr =  prop.ctr,
                          allocRatio = allocRatio,
                          CommonStdDev = CommonStdDev)
  SigmaZ_Stage <- matrix(outgetSigma$SigmaZ$EP1, nrow = 4)
  sigmaZ_Stage1 <- SigmaZ_Stage[c(1, 2), c(1, 2)]
  sigmaZ_Stage2 <- SigmaZ_Stage[c(3, 4), c(3, 4)]

  # getPlanCorrelation Results
  outgetPlanCorr <- getPlanCorrelation(
    nHypothesis = 2,
    EpType = EpType,
    SS_Incr = SS_Incr,
    Arms.std.dev = sigma,
    prop.ctr = prop.ctr,
    test.type = "Dunnett",
    CommonStdDev = CommonStdDev
  )

  corr_stage1 <- matrix(outgetPlanCorr$Stage1, nrow = 2)
  corr_stage2 <- matrix(outgetPlanCorr$Stage2, nrow = 2)

  expect_equal(object = corr_stage1, expected = sigmaZ_Stage1)
  expect_equal(object = corr_stage2, expected = sigmaZ_Stage2)
})
