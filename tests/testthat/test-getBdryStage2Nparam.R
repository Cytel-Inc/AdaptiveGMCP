

test_that('Test Non-Parametric Stage-2 plan boundary',{

  #Test Case-1: Spending Function: LD(OF) Validate with East Result
  ej <- 0.025
  aj1 <- 0.00152532274155114
  ss1 <- 50
  ss2 <- 100
  out <- getBdryStage2Nparam(ej=ej,aj1=aj1,ss1=ss1,ss2=ss2)
  benchmark <- 0.0244997717772641
  expect_equal(object = out, expected = benchmark, tolerance = 1E-9)

  #Test Case-2: extreme alpha Validate with East Result
  ej <- 0.001
  aj1 <- 3.26335376173607E-06
  ss1 <- 50
  ss2 <- 100
  out <- getBdryStage2Nparam(ej=ej,aj1=aj1,ss1=ss1,ss2=ss2)
  benchmark <- 0.000998794890168816
  expect_equal(object = out, expected = benchmark, tolerance = 1E-6)

  #Test Case-3: No Efficacy at stage-1
  ej <- 0.025
  aj1 <- 0
  ss1 <- 30
  ss2 <- 100
  out <- getBdryStage2Nparam(ej=ej,aj1=aj1,ss1=ss1,ss2=ss2)
  expect_equal(object = out, expected = ej, tolerance = 1E-9)

  #Test Case-4: significance level 0
  ej <- 0
  aj1 <- 0.00152532274155114
  ss1 <- 50
  ss2 <- 100
  out <- getBdryStage2Nparam(ej=ej,aj1=aj1,ss1=ss1,ss2=ss2)
  benchmark <- 0
  expect_equal(object = out, expected = benchmark, tolerance = 1E-9)

  #Test Case-4: Test Error if sample size inputs are wrong
  ej <- 0.025
  aj1 <- 0.00152532274155114
  expect_error(getBdryStage2Nparam(ej=ej,aj1=aj1,ss1=0,ss2=100))
  expect_error(getBdryStage2Nparam(ej=ej,aj1=aj1,ss1=100,ss2=0))

})
