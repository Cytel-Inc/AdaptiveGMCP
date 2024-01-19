



test_that('Test PCER computations',{
  #Test Case1: VConsistency with East 2-sample Muller-Schefer Benchmarks
  a2 <- 0.0244997717772641
  ss1 <- 50
  ss2 <- 100
  p1 <- 0.0454757
  benchmark <- 0.137062084206029
  out <- getPCER(a2=a2, p1=p1, ss1=ss1, ss2=ss2)
  expect_equal(object = out, expected = benchmark, tolerance = 1E-6)

  #Test Case2: Consistency with East 2-sample Muller-Schefer Benchmarks
  a2 <- 0.0244997717772641
  ss1 <- 110
  ss2 <- 219
  p1 <- 0.03963089
  benchmark <- 0.152051635710302
  out <- getPCER(a2=a2, p1=p1, ss1=ss1, ss2=ss2)
  expect_equal(object = out, expected = benchmark, tolerance = 1E-3)

  #Test Case3: consistency p1=0 => PCER = 1, p1=1 => PCER = 0
  a2 <- 0.0244997717772641
  ss1 <- 110
  ss2 <- 219
  p1 <- c(0,1)
  benchmark <- c(1,0)
  out <- getPCER(a2=a2, p1= p1, ss1=ss1, ss2=ss2)
  expect_equal(object = out, expected = benchmark)

  #Test Case4: expect error if ss1 = 0 or ss2 = 0
  expect_error(getPCER(a2=0.0244, p1= 0.03963089, ss1=0, ss2=100))
  expect_error(getPCER(a2=0.0244, p1= 0.03963089, ss1=50, ss2=0))


})
