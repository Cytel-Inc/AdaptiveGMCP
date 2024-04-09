test_that("Test Plan Parametric Boundary(CER Method)", {
  # Test Case1: Validate with East MAMS-GS boundary
  # Arms : 3, Looks: 2, Hypothesis: 2
  EpType <- list('EP1'='Continuous')
  alpha <- 0.025
  nLooks <- 2
  info_frac <- c(0.5, 1)
  wJ <- c(0.5, 0.5)
  sigma <- list("EP1" = c(1, 1, 1))
  prop.ctr <- NA
  allocRatio <- c(1, 1, 1)
  SS_Cum <- matrix(c(
    50, 50, 50,
    100, 100, 100
  ), nrow = 2, byrow = T)

  Sigma <- getSigma(EpType = EpType,
                    SS_Cum = SS_Cum,
                    sigma = sigma,
                    prop.ctr = prop.ctr,
                    allocRatio = allocRatio)
  gIDX <- 1
  hIDX <- c(1, 2)
  typeOfDesign <- "asOF"
  Scale <- "Score"
  out <- getPlanParmBdry(gIDX, hIDX, alpha, nLooks, info_frac, wJ, Sigma, typeOfDesign, Scale)
  pBdry <- rbind(out$Stage1Bdry, out$Stage2Bdry)

  # Benchmarks From East
  benchmark <- matrix(c(
    0.00078229234168911, 0.00078229234168911,
    0.0131202663476177, 0.0131202663476177
  ), nrow = 2, byrow = T)

  expect_equal(object = pBdry, expected = benchmark, tolerance = 1E-2)

  #--------------------------------------------------------------------------------
  # Test Case2: Validate with East MAMS-GS boundary unbalanced design
  # Arms : 3, Looks: 2, Hypothesis: 2
  EpType <- list('EP1'='Continuous')
  alpha <- 0.025
  nLooks <- 2
  info_frac <- c(0.5, 1)
  wJ <- c(0.5, 0.5)
  sigma <- list("EP1" = c(1.1, 1.2, 1.3))
  allocRatio <- c(1, 2, 1)
  SS_Cum <- matrix(c(
    50, 100, 50,
    100, 200, 100
  ), nrow = 2, byrow = T)

  Sigma <- getSigma(EpType = EpType,
                    SS_Cum = SS_Cum,
                    sigma = sigma,
                    allocRatio = allocRatio)
  gIDX <- 1
  hIDX <- c(1, 2)
  typeOfDesign <- "asOF"
  Scale <- "Score"
  out <- getPlanParmBdry(gIDX, hIDX, alpha, nLooks, info_frac, wJ, Sigma, typeOfDesign, Scale)
  pBdry <- rbind(out$Stage1Bdry, out$Stage2Bdry)

  # Benchmarks From East
  benchmark <- matrix(c(
    0.00141441817289648, 0.000126137793428335,
    0.0199601761219444, 0.00589281584231451
  ), nrow = 2, byrow = T)

  expect_equal(object = pBdry, expected = benchmark, tolerance = 1E-2)
  # Note: unequal sigma or allocation will not create different boundaries for different treatments with getPlanParmBdry

  #----------------------------------------------------------------------
  # Test Case3: Consistency with multiple endpoints
  # Arms : 3, Looks: 2, Endpoints: 3
  EpType <- list('EP1'='Continuous',
                 'EP2'='Continuous',
                 'EP3'='Continuous')
  alpha <- 0.025
  nLooks <- 2
  info_frac <- c(0.5, 1)
  wJ <- c(0.5, 0.5)
  sigma <- list(
    "EP1" = c(1, 1, 1),
    "EP2" = c(1, 1, 1),
    "EP3" = c(1, 1, 1)
  )
  allocRatio <- c(1, 1, 1)
  SS_Cum <- matrix(c(
    50, 50, 50,
    100, 100, 100
  ), nrow = 2, byrow = T)

  Sigma <- getSigma(EpType = EpType,
                    SS_Cum = SS_Cum,
                    sigma = sigma,
                    allocRatio = allocRatio)

  # EP-1
  gIDX <- 1
  hIDX <- c(1, 2)
  typeOfDesign <- "asOF"
  Scale <- "Score"
  out <- getPlanParmBdry(gIDX, hIDX, alpha, nLooks, info_frac, wJ, Sigma, typeOfDesign, Scale)
  pBdryEP1 <- rbind(out$Stage1Bdry, out$Stage2Bdry)

  # EP-2
  gIDX <- 2
  hIDX <- c(1, 2)
  typeOfDesign <- "asOF"
  Scale <- "Score"
  out <- getPlanParmBdry(gIDX, hIDX, alpha, nLooks, info_frac, wJ, Sigma, typeOfDesign, Scale)
  pBdryEP2 <- rbind(out$Stage1Bdry, out$Stage2Bdry)

  # EP-3
  gIDX <- 3
  hIDX <- c(1, 2)
  typeOfDesign <- "asOF"
  Scale <- "Score"
  out <- getPlanParmBdry(gIDX, hIDX, alpha, nLooks, info_frac, wJ, Sigma, typeOfDesign, Scale)
  pBdryEP3 <- rbind(out$Stage1Bdry, out$Stage2Bdry)

  expect_equal(
    object = all(pBdryEP1 == pBdryEP2) & all(pBdryEP1 == pBdryEP3),
    expected = TRUE
  )

  #-----------------------------------------------------------------------------
  # Test Case4: Validate with East MAMS-GS boundary
  # Arms : 3, Looks: 2, Hypothesis: 2
  EpType <- list('EP1'='Continuous')
  alpha <- 0.025
  nLooks <- 2
  info_frac <- c(0.5, 1)
  wJ <- c(0.5, 0.5)
  sigma <- list("EP1" = c(1, 1, 1))
  allocRatio <- c(1, 1, 1)
  SS_Cum <- matrix(c(
    50, 50, 50,
    100, 100, 100
  ), nrow = 2, byrow = T)

  Sigma <- getSigma(EpType = EpType,
                    SS_Cum = SS_Cum,
                    sigma = sigma,
                    prop.ctr =  prop.ctr,
                    allocRatio = allocRatio)
  gIDX <- 1
  hIDX <- c(1, 2)
  typeOfDesign <- "noEarlyEfficacy"
  Scale <- "Score"
  out <- getPlanParmBdry(gIDX = gIDX,
                         hIDX = hIDX,
                         alpha = alpha,
                         nLooks = nLooks,
                         info_frac = info_frac,
                         wJ = wJ,
                         Sigma =Sigma,
                         typeOfDesign = typeOfDesign,
                         Scale = Scale)
  pBdry <- rbind(out$Stage1Bdry, out$Stage2Bdry)

  # Benchmarks From East
  benchmark <- matrix(c(
    0, 0,
    0.0134479578529462, 0.0134479578529462
  ), nrow = 2, byrow = T)

  expect_equal(object = pBdry, expected = benchmark, tolerance = 1E-2)
})
