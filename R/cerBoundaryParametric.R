# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------- -
########### probability of crossing the boundary at stage-2##############
# cJ1    : The stage-1 critical point for the intersection hypothesis HJ
# wJ     : The weight for the intersection hypothesis HJ
# sigmaZ : The co-variance matrix of cummulative z-statistics (2 stages combined)
# underNull : TRUE if the probability is under null
# Returns: The Probability of rejecting in atleast one primary hypothesis at stage1
#------------------------------------------------------------------------- -
exitProbStage2 <- function(gIDX, hIDX, cJ2, cJ1, wJ, Sigma, Scale, 
          mvtnorm_algo, underNull = TRUE) {
  # Make Seed restricted to local##
  old <- .Random.seed
  on.exit({
    .Random.seed <<- old
  })
  set.seed(200295)
  #################################

  # browser()

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
    # stage-2 index
    sigmaIDX <- c(sIDX, sIDX + (ncol(sigmaZ) / 2))
    if (underNull) {
      mu_z <- rep(0, length(hIDX))
    }
    sigma <- sigmaZ[sigmaIDX, sigmaIDX]
    upper <- c(qnorm(1 - wJ[hIDX] * cJ1), qnorm(1 - wJ[hIDX] * cJ2))
    lower <- -Inf

    # Use dimension-based algorithm selected in simMAMSMEP()
    prob <- mvtnorm::pmvnorm(
      lower = lower, upper = upper, mean = mu_z, sigma = sigma, 
      algorithm = mvtnorm_algo
    )[1]
    (1 - prob) # Under null this should be cummulative alpha for that look
  } else if (Scale == "Score") {
    # stage-2 index
    sigmaIDX <- c(sIDX, sIDX + (ncol(sigmaS) / 2))

    if (underNull) {
      mu_s <- rep(0, length(hIDX))
    }
    sigma <- sigmaS[sigmaIDX, sigmaIDX]
    upper <- c(
      qnorm(1 - wJ[hIDX] * cJ1) * sqrt(infoMatrix[sIDX, 1]),
      qnorm(1 - wJ[hIDX] * cJ2) * sqrt(infoMatrix[sIDX, 2])
    )
    lower <- -Inf

    # Use dimension-based algorithm selected in simMAMSMEP()
    prob <- mvtnorm::pmvnorm(
      lower = lower, upper = upper, mean = mu_s, sigma = sigma, 
      algorithm = mvtnorm_algo
    )[1]
    (1 - prob) # Under null this should be cummulative alpha for that look
  }
}

#------------------------------------------------------------------------- -
#################### Conditional boundary crossing Prob ##########################
# cJ2    : The stage-2 critical point for the intersection hypothesis HJ
# p1     : stage-1 raw p-values for the subset Jh
# w      : The weight for the for the subset Jh
# sigmaSincr : The co-variance matrix for stage-2
# Conditional: TRUE = Conditional on the stage-1 p-values(p1)
#------------------------------------------------------------------------- -


exitProbStage2Cond <- function(cJ2, p1, w, InfoMatrix, stage2sigmaS, 
          mvtnorm_algo, Conditional = TRUE) {
  # Make Seed restricted to local##
  old <- .Random.seed
  on.exit({
    .Random.seed <<- old
  })
  set.seed(200295)
  #################################

  if (Conditional) {
    p1[p1==1] = 0.999999
    upper <- sqrt(InfoMatrix[, 2]) * qnorm(1 - w * cJ2) - sqrt(InfoMatrix[, 1]) * qnorm(1 - p1)
  } else {
    upper <- sqrt(InfoMatrix[, 2]) * qnorm(1 - w * cJ2)
  }
  lower <- -Inf

  # Use dimension-based algorithm selected in simMAMSMEP()
  1 - mvtnorm::pmvnorm(
    lower = lower, upper = upper, sigma = stage2sigmaS, algorithm = mvtnorm_algo
  )[1]
}
