# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

########### Boundary Computation Parametric Method ###########
#------------------------------------------------------------------------- -
########### probability of crossing the boundary at stage 1##############
# cJ1    : the stage-1 critical point for the intersection hypothesis HJ
# wJ     : the weight for the intersection hypothesis HJ
# sigmaZ : the co-variance matrix of cummulative z-statistics (2 stages combined)
# underNull : TRUE if the probability is under null
# Returns: The Probability of rejecting in atleast one primary hypothesis at stage1
#------------------------------------------------------------------------- -
exitProbParamStage1 <- function(gIDX, hIDX, cJ1, wJ, Sigma, Scale, mvtnorm_algo, underNull = TRUE) {
  # Make Seed restricted to local##
  old <- .Random.seed
  on.exit({
    .Random.seed <<- old
  })
  set.seed(200295)
  #################################

  SigmaZ <- Sigma$SigmaZ
  SigmaS <- Sigma$SigmaS
  InfoMatrix <- Sigma$InfoMatrix

  sigmaZ <- SigmaZ[[gIDX]]
  sigmaS <- SigmaS[[gIDX]]
  infoMatrix <- InfoMatrix[[gIDX]]

  nhyp <- length(wJ) # Number of hypothesis
  ngrps <- length(SigmaZ) # number of parametric groups

  sIDX <- rep(1:(nhyp / ngrps), ngrps)[hIDX] # stage-1 sigma index

  if (Scale == "Z") {
    if (underNull) {
      mu_z <- rep(0, length(sIDX))
    }
    sigma <- sigmaZ[sIDX, sIDX]
    upper <- qnorm(1 - wJ[hIDX] * cJ1)
    lower <- -Inf

    # Use dimension-based algorithm selected in simMAMSMEP()
    1 - mvtnorm::pmvnorm(
      lower = lower, upper = upper, mean = mu_z, sigma = sigma, algorithm = mvtnorm_algo
    )[1]
  } else if (Scale == "Score") {
    if (underNull) {
      mu_s <- rep(0, length(sIDX))
    }
    sigma <- sigmaS[sIDX, sIDX]
    upper <- qnorm(1 - wJ[hIDX] * cJ1) * sqrt(infoMatrix[sIDX, 1])
    lower <- -Inf

    # Use dimension-based algorithm selected in simMAMSMEP()
    1 - mvtnorm::pmvnorm(
      lower = lower, upper = upper, mean = mu_s, sigma = sigma, algorithm = mvtnorm_algo
    )[1]
  }
}

exitProbNParamStage1 <- function(cJ1, wJ, hIDX) {
  # we are taking some dummy value for now will change later once code works
  upper = qnorm(1 - wJ[hIDX] * cJ1)
  1 - pnorm(upper)
}
