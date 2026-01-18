# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

# Component-wise PCER
getPCER <- function(a2, p1, ss1, ss2, t1) {
  if (ss1 == 0 || ss2 == 0) stop("Error: ss1 or ss2 is zero | function:getPCER")
  # local Seed##
  set.seed(200295)
  #################################
  # r <- sqrt(ss1 / ss2)
  r <- sqrt(t1)

  if(p1 == 0){
    #if the stage-1 observed p-value is 0 then the PCER is 1
    1
  }else if(p1 == 1){
    #if the stage-1 observed p-value is 1 then the PCER is 0
    0
  }else{
    1 - pnorm(qnorm(1 - a2),
              mean = r * qnorm(1 - p1),
              sd = sqrt(1 - r^2))
  }
}

getPCER2 <- function(p1, ss1, ss2, cJ2, wJh, t1) {
  getPCER(a2 = cJ2 * wJh, p1 = p1, ss1 = ss1, ss2 = ss2, t1 = t1)
}

exitProbStage2Nparam2 <- function(cJ2, cJ1, ss1, ss2, wJ, hIDX, t1) # aj2: stage 2 bdry, aj1: stage-1 boundary
{
  aj2 = cJ2 * wJ[hIDX]
  aj1 = cJ1 * wJ[hIDX]
  # Set local seed##
  set.seed(200295)
  #################################
  upper <- c(qnorm(1 - aj1), Inf)
  lower <- c(-Inf, ifelse(1 - aj2 >= 0, qnorm(1 - aj2), -Inf))
  # Fixed a bug that was causing minor error in boundary computation 
  # due to the info fraction calculated as a ratio of sample sizes rather than
  # using the input info fraction directly.
  r <- sqrt(t1) # sqrt(ss1 / ss2)
  sigma <- matrix(c(1, r, r, 1), nrow = 2)
  prob <- mvtnorm::pmvnorm(lower = lower, upper = upper, sigma = sigma)[1]
  return(prob + aj1)
}
