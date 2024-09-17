########### Boundary Computation Parametric Method ###########
#------------------------------------------------------------------------- -
########### probability of crossing the boundary at stage 1##############
# cJ1    : the stage-1 critical point for the intersection hypothesis HJ
# wJ     : the weight for the intersection hypothesis HJ
# sigmaZ : the co-variance matrix of cummulative z-statistics (2 stages combined)
# underNull : TRUE if the probability is under null
# Returns: The Probability of rejecting in atleast one primary hypothesis at stage1
#------------------------------------------------------------------------- -
exitProbParamStage1 <- function(gIDX, hIDX, cJ1, wJ, Sigma, Scale, underNull = TRUE) {
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
    1 - mvtnorm::pmvnorm(
      lower = lower, upper = upper, mean = mu_z, sigma = sigma
    )[1]
  } else if (Scale == "Score") {
    if (underNull) {
      mu_s <- rep(0, length(sIDX))
    }
    sigma <- sigmaS[sIDX, sIDX]
    upper <- qnorm(1 - wJ[hIDX] * cJ1) * sqrt(infoMatrix[sIDX, 1])
    lower <- -Inf
    1 - mvtnorm::pmvnorm(
      lower = lower, upper = upper, mean = mu_s, sigma = sigma
    )[1]
  }
}

exitProbNParamStage1 <- function(cJ1, wJ, hIDX) {
  # we are taking some dummy value for now will change later once code works
  upper = qnorm(1 - wJ[hIDX] * cJ1)
  1 - pnorm(upper)
}

#------------------------------------------------------------------------- -
#################### Search for stage-1 boundary ##########################
# alpha1    : alpha spent at stage 1
# wJ     : The weight for the intersection hypothesis HJ
# sigmaZ : The co-variance matrix of cummulative z-statistics (2 stages combined)
# Returns: The critical point cJ1 for testing HJ at stage1
#------------------------------------------------------------------------- -
getBdryStage1 <- function(gIDX, hIDX, alpha1, wJ, Sigma, Scale) {
  minbdry <- 0#  0.00000000001
  maxbdry <- 1 / max(wJ[wJ != 0])
  bdry1 <- function(x) {
    extProbParam <- exitProbParamStage1(gIDX = gIDX, hIDX = hIDX, cJ1 = x, wJ = wJ, Sigma = Sigma, Scale = Scale, underNull = TRUE)
    # cat('cJ1 :',x,'extProb:',extProb,'\n')
    extProbParam - alpha1
  }
  uniroot(f = bdry1, interval = c(minbdry, maxbdry), tol = 1E-16)$root
}

#------------------------------------------------------------------------- -
########### probability of crossing the boundary at stage-2##############
# cJ1    : The stage-1 critical point for the intersection hypothesis HJ
# wJ     : The weight for the intersection hypothesis HJ
# sigmaZ : The co-variance matrix of cummulative z-statistics (2 stages combined)
# underNull : TRUE if the probability is under null
# Returns: The Probability of rejecting in atleast one primary hypothesis at stage1
#------------------------------------------------------------------------- -
exitProbParamStage2 <- function(gIDX, hIDX, cJ2, cJ1, wJ, Sigma, Scale, underNull = TRUE) {
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
    # stage-2 index
    sigmaIDX <- c(sIDX, sIDX + (ncol(sigmaZ) / 2))
    if (underNull) {
      mu_z <- rep(0, length(hIDX))
    }
    sigma <- sigmaZ[sigmaIDX, sigmaIDX]
    upper <- c(qnorm(1 - wJ[hIDX] * cJ1), qnorm(1 - wJ[hIDX] * cJ2))
    lower <- -Inf
    prob <- mvtnorm::pmvnorm(
      lower = lower,
      upper = upper,
      mean = mu_z,
      sigma = sigma,
      Seed = 200295
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
    prob <- mvtnorm::pmvnorm(
      lower = lower,
      upper = upper,
      mean = mu_s,
      sigma = sigma
    )[1]
    (1 - prob) # Under null this should be cummulative alpha for that look
  }
}

exitProbNParamStage2<- function(aj2, aj1, ss1, ss2) # aj2: stage 2 bdry, aj1: stage-1 boundary
  {
    # Set local seed##
    set.seed(200295)
    #################################
    upper <- c(qnorm(1 - aj1), Inf)
    lower <- c(-Inf, ifelse(1 - aj2 >= 0, qnorm(1 - aj2), -Inf))
    r <- sqrt(ss1 / ss2)
    sigma <- matrix(c(1, r, r, 1), nrow = 2)
    prob <- mvtnorm::pmvnorm(lower = lower, upper = upper, sigma = sigma)[1]
    return(prob + aj1)
  }


#------------------------------------------------------------------------- -
#################### Search for stage-2 boundary ##########################
# alpha  : planned alpha
# cJ1    : The stage-1 critical point for the intersection hypothesis HJ
# exitProb1 : Boundary crossing probability at stage-1
# wJ     : The weight for the intersection hypothesis HJ
# sigmaZ : The co-variance matrix of cummulative z-statistics (2 stages combined)
# Returns: The critical point cJ2 for testing HJ at stage2
#------------------------------------------------------------------------- -
getBdryStage2 <- function(gIDX, hIDX, alpha, cJ1, wJ, Sigma, Scale = Scale) {
  #browser()
  minbdry <-  0.00000000001
  maxbdry <- 1 / max(wJ[wJ != 0])
  bdry2 <- function(x) {
    extProbParam2 <- exitProbParamStage2(gIDX = gIDX, hIDX = hIDX, cJ2 = x, cJ1 = cJ1, wJ = wJ, Sigma = Sigma, Scale = Scale, underNull = TRUE)
    extProbParam2 - alpha
  }
  uniroot(f = bdry2, interval = c(minbdry, maxbdry), tol = 1E-16)$root
}


