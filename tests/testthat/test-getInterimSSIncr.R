

test_that("Test Per Look Simulation Sample Size(Incr.)",{
  #Arms: 4, Looks: 2
  #Test Case 1: No Implicit SSR
  lookID <- 2
  PlanSSIncr <- matrix(c(51,52,53,54,
                         51,52,53,54), nrow = 2, byrow = T)
  ArmsPresent <- c(T,T,T,F)
  ArmsRetained <- c(F,F,F,F)
  Arms.alloc.ratio <- c(1,1,1,1)
  ImplicitSSR <- 'None'

  Stage2SSIncr <- getInterimSSIncr(lookID = lookID, PlanSSIncr = PlanSSIncr,
                   ArmsPresent = ArmsPresent, ArmsRetained = ArmsRetained,
                   Arms.alloc.ratio = Arms.alloc.ratio, ImplicitSSR = ImplicitSSR)

  benchmark <- c(51,52,53,NA)
  expect_equal(object = Stage2SSIncr, expected = benchmark)
  expect_equal(sum(Stage2SSIncr, na.rm = T) == sum(PlanSSIncr[2,]), expected = F)
  #-------------------------------------------------------------------

  #Test Case 2: Implicit SSR only from the dropped arms due to treatment selection
  lookID <- 2
  PlanSSIncr <- matrix(c(51,52,53,54,
                         51,52,53,54), nrow = 2, byrow = T)
  ArmsPresent <- c(T,T,F,F)
  ArmsRetained <- c(F,F,F,T)
  Arms.alloc.ratio <- c(1,1,1,1)
  ImplicitSSR <- 'Selection'

  Stage2SSIncr <- getInterimSSIncr(lookID = lookID, PlanSSIncr = PlanSSIncr,
                                   ArmsPresent = ArmsPresent, ArmsRetained = ArmsRetained,
                                   Arms.alloc.ratio = Arms.alloc.ratio, ImplicitSSR = ImplicitSSR)

  benchmark <- c(51+54/2,52+54/2,NA,NA)
  expect_equal(object = Stage2SSIncr, expected = benchmark)
  expect_equal(sum(Stage2SSIncr, na.rm = T) == sum(PlanSSIncr[2,]), expected = F)

  #--------------------------------------------------------------------
  #Test Case 3: Implicit SSR from all dropped arms to continuing arms
  lookID <- 2
  PlanSSIncr <- matrix(c(51,52,56,54,
                         51,52,56,54), nrow = 2, byrow = T)
  ArmsPresent <- c(T,T,F,F)
  ArmsRetained <- c(F,F,F,T)
  Arms.alloc.ratio <- c(1,1,1,1)
  ImplicitSSR <- 'All'

  Stage2SSIncr <- getInterimSSIncr(lookID = lookID, PlanSSIncr = PlanSSIncr,
                                   ArmsPresent = ArmsPresent, ArmsRetained = ArmsRetained,
                                   Arms.alloc.ratio = Arms.alloc.ratio, ImplicitSSR = ImplicitSSR)

  benchmark <- c(51+(56+54)/2,52+(56+54)/2,NA,NA)
  expect_equal(object = Stage2SSIncr, expected = benchmark)
  expect_equal(sum(Stage2SSIncr, na.rm = T) == sum(PlanSSIncr[2,]), expected = T)

  #------------------------------------------------------------------
  #Test Case 4: Expect Error
  lookID <- 2
  PlanSSIncr <- matrix(c(51,52,56,54,
                         51,52,56,54), nrow = 2, byrow = T)
  ArmsPresent <- c(T,T,F)
  ArmsRetained <- c(F,F,F)
  Arms.alloc.ratio <- c(1,1,1)
  ImplicitSSR <- 'All'

  expect_error(getInterimSSIncr(lookID = lookID, PlanSSIncr = PlanSSIncr,
                                ArmsPresent = ArmsPresent, ArmsRetained = ArmsRetained,
                                Arms.alloc.ratio = Arms.alloc.ratio, ImplicitSSR = ImplicitSSR))

  lookID <- 2
  PlanSSIncr <- matrix(c(51,52,56,
                         51,52,56), nrow = 2, byrow = T)
  ArmsPresent <- c(F,F,F)
  ArmsRetained <- c(T,T,T)
  Arms.alloc.ratio <- c(1,1,1)
  ImplicitSSR <- 'All'

  expect_error(getInterimSSIncr(lookID = lookID, PlanSSIncr = PlanSSIncr,
                                ArmsPresent = ArmsPresent, ArmsRetained = ArmsRetained,
                                Arms.alloc.ratio = Arms.alloc.ratio, ImplicitSSR = ImplicitSSR))


})
