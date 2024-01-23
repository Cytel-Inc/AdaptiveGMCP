


test_that("Test Computations of per-look Summary Statistics Computations",{

  #Test Case1: Verify Look1 summary statistics with Excel Benchmarks(from the subject data)
  SimSeed <- 100
  simID <- 1
  lookID <- 1
  Arms.Mean <- list('EP1'=c(0,0.1,0.4))
  Arms.std.dev <- list('EP1'=c(1,1,1))
  Arms.alloc.ratio <- c(1,1,1)
  Arms.SS <- c(56,56,56)
  ArmsPresent <- c(T,T,T)
  HypoPresent <- c(T,T)
  HypoMap <- data.frame('Hypothesis'=c('H1','H2'),
                        'Groups'=1,
                        'Control'=1,
                        'Treatment'=c(2,3))
  Cumulative <- F

  Stage1Response <- genIncrLookSummaryDOM(SimSeed=SimSeed,
                        simID=simID,
                        lookID=lookID,
                        Arms.Mean=Arms.Mean,
                        Arms.std.dev=Arms.std.dev,
                        Arms.alloc.ratio=Arms.alloc.ratio,
                        Arms.SS=Arms.SS,
                        ArmsPresent=ArmsPresent,
                        HypoPresent=HypoPresent,
                        HypoMap =  HypoMap)

  SummStat1 <- getPerLookTestStatDOM(simID = simID,
                                    lookID = lookID,
                                    IncrLookSummaryDOM=Stage1Response,
                                    Cumulative = Cumulative)

  delta_benchmark <- c(-0.14570479,	0.175890224)
  SE_benchmark <- c(0.189234032,	0.191108795)
  TestStat_benchmark <- c(-0.769971389,	0.920366978)
  pValue_benchmark <- 1-pnorm(TestStat_benchmark)

  delta <- unlist(SummStat1[,grep('Delta', names(SummStat1))])
  SE <- unlist(SummStat1[,grep('StdError', names(SummStat1))])
  TestStat <- unlist(SummStat1[,grep('TestStat', names(SummStat1))])
  pValue <- unlist(SummStat1[,grep('RawPvalues', names(SummStat1))])
  names(delta) <- names(SE) <- names(TestStat) <- names(pValue) <- NULL

  expect_equal(object = delta, expected = delta_benchmark)
  expect_equal(object = SE, expected = SE_benchmark)
  expect_equal(object = TestStat, expected = TestStat_benchmark)
  expect_equal(object = pValue, expected = pValue_benchmark)
  #-----------------------------------------------------------------------------------------------

  #Test Case2: Verify Look2 summary statistics(Incr.) with Excel Benchmarks(from the subject data)
  lookID <- 2
  Arms.SS <- c(56,56,56)
  ArmsPresent <- c(T,T,T)
  HypoPresent <- c(T,T)
  Cumulative <- F
  Stage2Response <- genIncrLookSummaryDOM(SimSeed=SimSeed,
                                          simID=simID,
                                          lookID=lookID,
                                          Arms.Mean=Arms.Mean,
                                          Arms.std.dev=Arms.std.dev,
                                          Arms.alloc.ratio=Arms.alloc.ratio,
                                          Arms.SS=Arms.SS,
                                          ArmsPresent=ArmsPresent,
                                          HypoPresent=HypoPresent,
                                          HypoMap =  HypoMap)

  SummStat2Incr <- getPerLookTestStatDOM(simID = simID,
                                    lookID = lookID,
                                    IncrLookSummaryDOM=Stage2Response,
                                    Cumulative = Cumulative)

  delta_benchmark <- c(0.187688947,	0.454780878)
  SE_benchmark <- c(0.167115243,	0.171058961)
  TestStat_benchmark <- c(1.123110874,	2.658620608)
  pValue_benchmark <- 1-pnorm(TestStat_benchmark)

  delta <- unlist(SummStat2Incr[,grep('Delta', names(SummStat2Incr))])
  SE <- unlist(SummStat2Incr[,grep('StdError', names(SummStat2Incr))])
  TestStat <- unlist(SummStat2Incr[,grep('TestStat', names(SummStat2Incr))])
  pValue <- unlist(SummStat2Incr[,grep('RawPvalues', names(SummStat2Incr))])
  names(delta) <- names(SE) <- names(TestStat) <- names(pValue) <- NULL

  expect_equal(object = delta, expected = delta_benchmark)
  expect_equal(object = SE, expected = SE_benchmark)
  expect_equal(object = TestStat, expected = TestStat_benchmark)
  expect_equal(object = pValue, expected = pValue_benchmark)

  #-----------------------------------------------------------------------------

  #Test Case3: Verify Look2 summary statistics(Cum.) with Excel Benchmarks(from the subject data)
  Cumulative <- T
  SummStat2Cum <- getPerLookTestStatDOM(simID = simID,
                                         lookID = lookID,
                                         IncrLookSummaryDOM = Stage2Response,
                                         IncrLookSummaryDOMPrev = Stage1Response,
                                         Cumulative = Cumulative)

  delta_benchmark <- c(0.020992078,	0.315335551)
  SE_benchmark <- c(0.126197841,	0.128012303)
  TestStat_benchmark <- c(0.166342612,	2.463322229)
  pValue_benchmark <- 1-pnorm(TestStat_benchmark)

  delta <- unlist(SummStat2Cum[,grep('Delta', names(SummStat2Cum))])
  SE <- unlist(SummStat2Cum[,grep('StdError', names(SummStat2Cum))])
  TestStat <- unlist(SummStat2Cum[,grep('TestStat', names(SummStat2Cum))])
  pValue <- unlist(SummStat2Cum[,grep('RawPvalues', names(SummStat2Cum))])
  names(delta) <- names(SE) <- names(TestStat) <- names(pValue) <- NULL

  expect_equal(object = delta, expected = delta_benchmark)
  expect_equal(object = SE, expected = SE_benchmark)
  expect_equal(object = TestStat, expected = TestStat_benchmark)
  expect_equal(object = pValue, expected = pValue_benchmark)
  #------------------------------------------------------------------------------

})
