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
exitProbStage1 <- function(gIDX, hIDX, cJ1, wJ, Sigma, Scale, underNull = TRUE) {
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
      lower = lower, upper = upper, mean = mu_z, sigma = sigma,
      algorithm = mvtnorm::Miwa(
        steps = 128,
        checkCorr = F,
        maxval = 1e3)
    )[1]
  } else if (Scale == "Score") {
    if (underNull) {
      mu_s <- rep(0, length(sIDX))
    }
    sigma <- sigmaS[sIDX, sIDX]
    upper <- qnorm(1 - wJ[hIDX] * cJ1) * sqrt(infoMatrix[sIDX, 1])
    lower <- -Inf
    1 - mvtnorm::pmvnorm(
      lower = lower, upper = upper, mean = mu_s, sigma = sigma,
      algorithm = Miwa(
        steps = 128,
        checkCorr = F,
        maxval = 1e3)
    )[1]
  }
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
    extProb <- exitProbStage1(gIDX = gIDX, hIDX = hIDX, cJ1 = x, wJ = wJ, Sigma = Sigma, Scale = Scale, underNull = TRUE)
    # cat('cJ1 :',x,'extProb:',extProb,'\n')
    extProb - alpha1
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
exitProbStage2 <- function(gIDX, hIDX, cJ2, cJ1, wJ, Sigma, Scale, underNull = TRUE) {
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
      algorithm = mvtnorm::Miwa(
        steps = 128,
        checkCorr = F,
        maxval = 1e3)
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
      sigma = sigma,
      algorithm = mvtnorm::Miwa(
        steps = 128,
        checkCorr = F,
        maxval = 1e3)
    )[1]
    (1 - prob) # Under null this should be cummulative alpha for that look
  }
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
  minbdry <-  0.00000000001
  maxbdry <- 1 / max(wJ[wJ != 0])
  bdry2 <- function(x) {
    extProb <- exitProbStage2(gIDX = gIDX, hIDX = hIDX, cJ2 = x, cJ1 = cJ1, wJ = wJ, Sigma = Sigma, Scale = Scale, underNull = TRUE)
    # cat('cJ2 :',x,'extProb:',extProb,'\n')
    extProb - alpha
  }
  uniroot(f = bdry2, interval = c(minbdry, maxbdry), tol = 1E-16)$root
}

#------------------------------------------------------------------------- -
#################### stage-1 & stage-2 boundary ##########################
# alpha  : planned alpha
# cJ1    : The stage-1 critical point for the intersection hypothesis HJ
# exitProb1 : Boundary crossing probability at stage-1
# wJ     : The weight for the intersection hypothesis HJ
# sigmaZ : The co-variance matrix of cummulative z-statistics (2 stages combined)
# Returns: The critical point cJ2 for testing HJ at stage2
#------------------------------------------------------------------------- -
getPlanParmBdry <- function(gIDX,
                            hIDX,
                            alpha,
                            nLooks,
                            info_frac,
                            wJ,
                            Sigma,
                            typeOfDesign,
                            deltaWT,
                            deltaPT1,
                            gammaA,
                            userAlphaSpending,
                            Scale) {

  if(typeOfDesign == "WT"){
    alpha1 <- rpact::getDesignGroupSequential(
      kMax = nLooks, alpha = alpha,
      informationRates = info_frac,
      typeOfDesign = typeOfDesign,
      deltaWT = deltaWT
    )$alphaSpent[1]
  }else if(typeOfDesign == "PT"){
    alpha1 <- rpact::getDesignGroupSequential(
      kMax = nLooks, alpha = alpha,
      informationRates = info_frac,
      typeOfDesign = typeOfDesign,
      deltaPT1 = deltaPT1
    )$alphaSpent[1]
  }else if(typeOfDesign == "asHSD" || typeOfDesign == "asKD"){
    alpha1 <- rpact::getDesignGroupSequential(
      kMax = nLooks, alpha = alpha,
      informationRates = info_frac,
      typeOfDesign = typeOfDesign,
      gammaA = gammaA
    )$alphaSpent[1]
  }else{
    alpha1 <- rpact::getDesignGroupSequential(
      kMax = nLooks, alpha = alpha,
      informationRates = info_frac,
      typeOfDesign = typeOfDesign
    )$alphaSpent[1]
  }


  if (nLooks == 1) {
    # Stage-1 boundary
    cJ1 <- getBdryStage1(gIDX = gIDX, hIDX = hIDX, alpha1 = alpha1, wJ = wJ, Sigma = Sigma, Scale = Scale)
    list("Stage1Bdry" = cJ1 * wJ, "Stage2Bdry" = NA)
  } else if (nLooks == 2) {
    # Stage-1 boundary
    cJ1 <- getBdryStage1(gIDX = gIDX, hIDX = hIDX, alpha1 = alpha1, wJ = wJ, Sigma = Sigma, Scale = Scale)
    # Stage-2 boundary
    cJ2 <- getBdryStage2(gIDX = gIDX, hIDX = hIDX, alpha = alpha, cJ1 = cJ1, wJ = wJ, Sigma = Sigma, Scale = Scale)
    list("Stage1Bdry" = cJ1 * wJ, "Stage2Bdry" = cJ2 * wJ)
  } else {
    print("Error in getPlanParmBdry: Boundary Computation is not available for Stages > 2")
  }
}

#------------------------------------------------------------------------- -


#------------------------------------------------------------------------- -
#################### Conditional boundary crossing Prob ##########################
# cJ2    : The stage-2 critical point for the intersection hypothesis HJ
# p1     : stage-1 raw p-values for the subset Jh
# w      : The weight for the for the subset Jh
# sigmaSincr : The co-variance matrix for stage-2
# Conditional: TRUE = Conditional on the stage-1 p-values(p1)
#------------------------------------------------------------------------- -


exitProbStage2Cond <- function(cJ2, p1, w, InfoMatrix, stage2sigmaS, Conditional = TRUE) {
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
  1 - mvtnorm::pmvnorm(
    lower = lower, upper = upper, sigma = stage2sigmaS
  )[1]
}


#------------------------------------------------------------------------- -
#################### Search for stage-1 boundary ##########################
# alpha1    : alpha spent at stage 1
# wJ     : The weight for the intersection hypothesis HJ
# sigmaZ : The co-variance matrix of cummulative z-statistics (2 stages combined)
# Returns: The critical point cJ1 for testing HJ at stage1
#------------------------------------------------------------------------- -
getStage2CondParamBdry <- function(cer, p1, w, InfoMatrix, stage2sigmaS, Conditional) {
  minbdry <- 0.0000001
  maxbdry <- 1 / max(w[w != 0])
  bdryCond <- function(x) {
    exitProb <- exitProbStage2Cond(
      cJ2 = x, p1 = p1, w = w, InfoMatrix = InfoMatrix,
      stage2sigmaS = stage2sigmaS, Conditional = Conditional
    )
    exitProb - cer
  }
  uniroot(f = bdryCond, interval = c(minbdry, maxbdry), tol = 1E-16)$root
}
