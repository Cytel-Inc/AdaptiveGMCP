# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

################# Computation of Stage-Wise covariance Matrix ###########

#------------------------------------------------------------------------ -
########### Elements of the Z-scale Covariance Matrix##############
# (i1,i2) : index for the hypothesis
# (k1,k2) : index for the looks
# sigma_0 : Standard deviation for the control arm
# sigma_trt : Standard deviation for the treatment arms
# ctrSS : stage-wise cumulative samples in control arm
# trtSS : stage-wise cumulative samples in treatment arms
# InfoMatrix : Fisher information matrix
# Returns: The Probability of rejecting in atleast one primary hypothesis at stage1
#------------------------------------------------------------------------- -
varCovZ <- function(EpType, i1, k1, i2, k2, sigma_0, sigma_trt, ctrProp, ctrSS, trtSS, InfoMatrix) {
  if (EpType == "Binary") {
    sigma_0 <- sigma_trt <- sqrt(ctrProp * (1 - ctrProp))
    if (i1 == i2 & k1 == k2) { # Variance Case-1
      1
    } else if (k1 == k2 & i1 != i2) # Covariance Case-2
      {
        sqrt(InfoMatrix[i1, k1] * InfoMatrix[i2, k2]) * (sigma_0^2 / ctrSS[k2])
      } else if (i1 == i2 & k1 != k2) # Covariance Case-3
      {
        sqrt(InfoMatrix[i1, k1] * InfoMatrix[i2, k2]) *
          ((sigma_trt^2 / trtSS[i2, max(k1, k2)]) + (sigma_0^2 / ctrSS[max(k1, k2)]))
      } else if (i1 != i2 || k1 != k2) # Covariance Case-4
      {
        sqrt(InfoMatrix[i1, k1] * InfoMatrix[i2, k2]) * (sigma_0^2 / ctrSS[max(k1, k2)])
      } else {
      "Error in covZ"
    }
  } else if (EpType == "Continuous") {
    if (i1 == i2 & k1 == k2) { # Variance Case-1
      1
    } else if (k1 == k2 & i1 != i2) # Covariance Case-2
      {
        sqrt(InfoMatrix[i1, k1] * InfoMatrix[i2, k2]) * (sigma_0^2 / ctrSS[k2])
      } else if (i1 == i2 & k1 != k2) # Covariance Case-3
      {
        sqrt(InfoMatrix[i1, k1] * InfoMatrix[i2, k2]) *
          ((sigma_trt[i2]^2 / trtSS[i2, max(k1, k2)]) + (sigma_0^2 / ctrSS[max(k1, k2)]))
      } else if (i1 != i2 || k1 != k2) # Covariance Case-4
      {
        sqrt(InfoMatrix[i1, k1] * InfoMatrix[i2, k2]) * (sigma_0^2 / ctrSS[max(k1, k2)])
      } else {
      "Error in covZ"
    }
  }
}

##### Computation of Covariance Matrix #############
# EpType : endpoint type
# SS_Cum : look-wise arm-wise cumulative samples
# sigma : list containing Arm-wise planned standard deviation for Continuous endpoints
# prop.ctr : list containing control proportion for Binary endpoints
# allocRatio : allocation ratio for the arms
#------------------------------------------------ -
getSigma <- function(EpType, SS_Cum, sigma, prop.ctr, allocRatio, CommonStdDev) {

  #Treatment standard Deviation is same is control
  #CommonStdDev flag is global variable set from simMAMSMEP(.), adaptGMCP_CER(.) function call
  if(CommonStdDev == T){
    for(sigIDX in 1:length(sigma))
      if(all(!is.na(sigma[[sigIDX]]))){
        sigma[[sigIDX]] <- sapply(sigma[[sigIDX]], function(x)sigma[[sigIDX]][1])
      }
  }

  ctrSS <- SS_Cum[, 1]
  trtSS <- t(SS_Cum[, -1]) # Column represents looks

  SigmaZ <- list() # Z-Scale sigma
  SigmaS <- list() # Score Scale sigma
  InfoMat <- list() # Information Matrix

  nEps <- length(EpType)


  for (i in 1:nEps) {
    if (EpType[[i]] == "Continuous") {
      nHypothesisEp <- length(sigma[[1]]) - 1 # Two equal dimension sigma matrix for two endpoints
      nLooksEp <- length(ctrSS)

      epSig <- sigma[[i]]
      sigma_0 <- epSig[1]
      sigma_trt <- epSig[-1]

      capLambda <- (sigma_0^2 + sigma_trt^2 / allocRatio[-1])^(-1)
    } else if (EpType[[i]] == "Binary") {
      nHypothesisEp <- length(allocRatio) - 1 # Two equal dimension sigma matrix for two endpoints
      nLooksEp <- length(ctrSS)
      pi_c <- prop.ctr[[i]]
      capLambda <- (1 / (pi_c * (1 - pi_c))) * (allocRatio[-1] / (1 + allocRatio[-1]))
    }


    InfoMatrix <- sapply(ctrSS, function(x) {
      x * capLambda
    }) # row=hypothesis, col=looks
    InfoMat[[paste("EP", i, sep = "")]] <- InfoMatrix

    ############### Computation of Z scale Sigma Matrix ################
    sigmaZ <- matrix(NA, nrow = nHypothesisEp * nLooksEp, ncol = nHypothesisEp * nLooksEp)
    hIDX <- rep(1:nHypothesisEp, nLooksEp)
    lIDX <- rep(1:nLooksEp, each = nHypothesisEp)

    for (l in 1:length(hIDX))
    {
      for (m in l:length(hIDX)) {
        if (EpType[[i]] == "Continuous") {
          sigmaZ[l, m] <- sigmaZ[m, l] <- varCovZ(
            EpType = "Continuous",
            i1 = hIDX[l], k1 = lIDX[l],
            i2 = hIDX[m], k2 = lIDX[m],
            sigma_0 = sigma_0, sigma_trt = sigma_trt,
            ctrSS = ctrSS, trtSS = trtSS, InfoMatrix = InfoMatrix
          )
        } else if (EpType[[i]] == "Binary") {
          sigmaZ[l, m] <- sigmaZ[m, l] <- varCovZ(
            EpType = "Binary",
            i1 = hIDX[l], k1 = lIDX[l],
            i2 = hIDX[m], k2 = lIDX[m],
            ctrProp = prop.ctr[[i]],
            ctrSS = ctrSS, trtSS = trtSS, InfoMatrix = InfoMatrix
          )
        }
      }
    }
    rownames(sigmaZ) <- colnames(sigmaZ) <- paste("Z", hIDX, lIDX, sep = "")
    SigmaZ[[paste("EP", i, sep = "")]] <- sigmaZ

    #----------------------------End of SigmaZ-------------------------

    ################## Computation of Score scale Sigma Matrix #########
    if (nLooksEp == 2) {
      l <- c(sqrt(InfoMatrix[, 1]), sqrt(InfoMatrix[, 2])) # this code is only for two stages
    } else if (nLooksEp == 1) {
      l <- c(sqrt(InfoMatrix[, 1])) # this code is only for  FSD
    } else {
      print("Error in getSigma: The conversion of Z-Score is not available for Stages > 2")
    }

    sigmaS <- matrix(NA, nrow = nrow(sigmaZ), ncol = ncol(sigmaZ))
    for (m in 1:nrow(sigmaS)) {
      for (n in m:ncol(sigmaS)) {
        sigmaS[m, n] <- sigmaS[n, m] <- l[m] * l[n] * sigmaZ[m, n]
      }
    }
    rownames(sigmaS) <- colnames(sigmaS) <- paste("S", hIDX, lIDX, sep = "")
    SigmaS[[paste("EP", i, sep = "")]] <- sigmaS
    #----------------------------End of SigmaS-------------------------
  }
  list("InfoMatrix" = InfoMat, "SigmaZ" = SigmaZ, "SigmaS" = SigmaS)
}



##### Covariance Matrix for CER & Stage-2 boundary computations#############
getStage2Sigma <- function(nHypothesis, EpType, nLooks, Sigma,
                           AllocSampleSize, allocRatio, sigma, prop.ctr,
                           Stage2AllocSampleSize, Stage2allocRatio, Stage2sigma,CommonStdDev) {
  #Treatment standard Deviation is same is control
  #CommonStdDev flag is global variable set from simMAMSMEP(.), adaptGMCP_CER(.) function call
  if(CommonStdDev == T){
    for(sigIDX in 1:length(sigma))
      if(all(!is.na(sigma[[sigIDX]]))){
        sigma[[sigIDX]] <- sapply(sigma[[sigIDX]], function(x)sigma[[sigIDX]][1])
      }
  }
  Stage2sigma <- sigma
  ############################################

  SigmaZIncr <- list() # Z-Scale Incremental sigma
  SigmaSIncr <- list() # Score Incremental Scale sigma
  Stage2SigmaZ <- list() # Z-Scale Modified sigma
  Stage2SigmaS <- list() # Score-Scale Modified sigma
  Stage2InfoMat <- list() # Stage-2 information matrix
  Stage2InfoMat_Incr <- list() # Stage-2 information matrix (incremental)

  ############# With planned Samples for CER computations##############
  nGrps <- length(Sigma$SigmaZ)
  for (i in 1:nGrps) # For all groups
  {
    # m <- ncol(Sigma$SigmaZ[[i]])
    # A <- getAmatrix(nrow = m/2, ncol = m)
    # SigmaZIncr[[names(Sigma$SigmaZ)[i]]] <- A %*% Sigma$SigmaZ[[i]] %*% t(A)
    # SigmaSIncr[[names(Sigma$SigmaS)[i]]] <- A %*% Sigma$SigmaS[[i]] %*% t(A)

    if (EpType[[i]] == "Continuous") {
      sigma_0 <- sigma[[i]][1]
      sigma_trt <- sigma[[i]][-1]
      capLambda <- (sigma_0^2 + sigma_trt^2 / allocRatio[-1])^(-1)
    } else if (EpType[[i]] == "Binary") {
      pi_c <- prop.ctr[[i]]
      capLambda <- (1 / (pi_c * (1 - pi_c))) * (allocRatio[-1] / (1 + allocRatio[-1]))
    }

    # AllocSampleSize : the planned sample size
    SSIncr <- as.numeric(AllocSampleSize[2, ]) - as.numeric(AllocSampleSize[1, ])
    ctrSS <- SSIncr[1]
    trtSS <- SSIncr[-1]
    InfoMatrix <- matrix(ctrSS * capLambda, ncol = 1)

    ########## Computation of  Z-scale Covariance Matrix #################
    k <- ncol(AllocSampleSize) - 1
    sigmaZ <- matrix(NA, nrow = k, ncol = k)

    for (l in 1:k)
    {
      for (m in l:k) {
        if (EpType[[i]] == "Continuous") {
          sigmaZ[l, m] <- sigmaZ[m, l] <- varCovZ(
            EpType = "Continuous",
            i1 = l, k1 = 1, i2 = m, k2 = 1,
            sigma_0 = sigma_0, sigma_trt = sigma_trt,
            ctrSS = ctrSS, trtSS = trtSS, InfoMatrix = InfoMatrix
          )
        } else if (EpType[[i]] == "Binary") {
          sigmaZ[l, m] <- sigmaZ[m, l] <- varCovZ(
            EpType = "Binary",
            i1 = l, k1 = 1, i2 = m, k2 = 1,
            ctrProp = pi_c,
            ctrSS = ctrSS, trtSS = trtSS, InfoMatrix = InfoMatrix
          )
        }
      }
    }
    SigmaZIncr[[names(Sigma$SigmaZ)[i]]] <- sigmaZ

    l <- c(sqrt(InfoMatrix[, 1])) # this code is only for  2-stage
    sigmaS <- matrix(NA, nrow = nrow(sigmaZ), ncol = ncol(sigmaZ))
    for (m in 1:nrow(sigmaS)) {
      for (n in m:ncol(sigmaS)) {
        sigmaS[m, n] <- sigmaS[n, m] <- l[m] * l[n] * sigmaZ[m, n]
      }
    }
    SigmaSIncr[[names(Sigma$SigmaZ)[i]]] <- sigmaS
  }

  ################ With modified Samples for Stage-2 boundary computations ######
  for (epIDX in 1:nGrps)
  {
    # The Cumulative InfoMatrix is needed for transformation(z to score)
    Stage2SSCum <- (Stage2AllocSampleSize[2, ])
    for (s in 1:length(Stage2SSCum)) {
      if (is.na(Stage2SSCum[s])) Stage2SSCum[s] <- AllocSampleSize[1, s]
    }
    Stage2allocRatio <- as.numeric(Stage2SSCum) / as.numeric(Stage2SSCum[1])

    if (EpType[[epIDX]] == "Continuous") {
      Stage2sigma_0 <- Stage2sigma[[epIDX]][1]
      Stage2sigma_trt <- Stage2sigma[[epIDX]][-1]
      Stage2capLambda <- (Stage2sigma_0^2 + Stage2sigma_trt^2 / Stage2allocRatio[-1])^(-1)
    } else if (EpType[[epIDX]] == "Binary") {
      pi_c <- prop.ctr[[epIDX]]
      Stage2capLambda <- (1 / (pi_c * (1 - pi_c))) * (Stage2allocRatio[-1] / (1 + Stage2allocRatio[-1]))
    }

    # The following adjustment is due to change in distribution(as the sample size modified)
    Stage2SSIncr <- as.numeric(Stage2AllocSampleSize[2, ]) - as.numeric(Stage2AllocSampleSize[1, ])
    Stage2ctrSS <- Stage2SSIncr[1]
    Stage2trtSS <- Stage2SSIncr[-1]
    Stage2InfoMatrix <- matrix(Stage2ctrSS * Stage2capLambda, ncol = 1)
    Stage2InfoMatrixCum <- matrix(as.numeric(Stage2SSCum[, 1]) * Stage2capLambda, ncol = 1)


    ########## Computation of  Z-scale Covariance Matrix #################
    k <- ncol(AllocSampleSize) - 1
    Stage2sigmaZ <- matrix(NA, nrow = k, ncol = k)

    for (l in 1:k)
    {
      for (m in l:k) {
        if (EpType[[epIDX]] == "Continuous") {
          Stage2sigmaZ[l, m] <- Stage2sigmaZ[m, l] <- varCovZ(
            EpType = "Continuous",
            i1 = l, k1 = 1, i2 = m, k2 = 1,
            sigma_0 = Stage2sigma_0, sigma_trt = Stage2sigma_trt,
            ctrSS = Stage2ctrSS, trtSS = Stage2trtSS,
            InfoMatrix = Stage2InfoMatrix
          )
        } else if (EpType[[epIDX]] == "Binary") {
          Stage2sigmaZ[l, m] <- Stage2sigmaZ[m, l] <- varCovZ(
            EpType = "Binary",
            i1 = l, k1 = 1, i2 = m, k2 = 1,
            ctrProp = pi_c,
            ctrSS = Stage2ctrSS, trtSS = Stage2trtSS,
            InfoMatrix = Stage2InfoMatrix
          )
        }
      }
    }
    Stage2SigmaZ[[names(Sigma$SigmaZ)[epIDX]]] <- Stage2sigmaZ
    Stage2InfoMat[[names(Sigma$SigmaZ)[epIDX]]] <- Stage2InfoMatrixCum
    Stage2InfoMat_Incr[[names(Sigma$SigmaZ)[epIDX]]] <- Stage2InfoMatrix

    ########## Computation of  Score-scale Covariance Matrix #################
    l <- c(sqrt(Stage2InfoMatrix[, 1])) # this code is only for  2-stage
    Stage2sigmaS <- matrix(NA, nrow = nrow(Stage2sigmaZ), ncol = ncol(Stage2sigmaZ))
    for (m in 1:nrow(Stage2sigmaS)) {
      for (n in m:ncol(Stage2sigmaS)) {
        Stage2sigmaS[m, n] <- Stage2sigmaS[n, m] <- l[m] * l[n] * Stage2sigmaZ[m, n]
      }
    }
    Stage2SigmaS[[names(Sigma$SigmaS)[epIDX]]] <- Stage2sigmaS
  }
  list(
    "SigmaZIncr" = SigmaZIncr,
    "SigmaSIncr" = SigmaSIncr,
    "Stage2SigmaZ" = Stage2SigmaZ,
    "Stage2SigmaS" = Stage2SigmaS,
    "Stage2InfoMatrixCum" = Stage2InfoMat,
    "Stage2InfoMatrixIncr" = Stage2InfoMat_Incr
  )
}



######### Difference Operator For Incremental Sigma############
getAmatrix <- function(nrow, ncol) {
  mat <- matrix(0, nrow = nrow, ncol = ncol)
  for (i in 1:nrow) {
    mat[i, i] <- -1
    mat[i, (i + ncol / 2)] <- 1
  }
  mat
}


######## Correlation for Combining p-value dunnett test#############
getPlanCorrelation <- function(nHypothesis, EpType, SS_Incr, Arms.std.dev, prop.ctr, test.type, CommonStdDev) {
  #Treatment standard Deviation is same is control
  #CommonStdDev flag is global variable set from simMAMSMEP(.), adaptGMCP_CER(.) function call
  if(CommonStdDev == T){
    for(sigIDX in 1:length(Arms.std.dev))
      if(all(!is.na(Arms.std.dev[[sigIDX]]))){
        Arms.std.dev[[sigIDX]] <- sapply(Arms.std.dev[[sigIDX]], function(x)Arms.std.dev[[sigIDX]][1])
      }
  }
  ############################################

  nEps <- length(EpType)
  nLooks <- nrow(SS_Incr)

  Sigma <- list()
  for (lkIDX in 1:nLooks) {
    SigmaZ <- list() # ith look Z-Scale sigma
    SS_lk <- as.numeric(SS_Incr[lkIDX, ])
    ctrSS <- SS_lk[1]
    trtSS <- SS_lk[-1] # Column represents looks
    allocRatio <- c(1, trtSS / ctrSS)

    #Expected Output : SS_lk(the correlation between test statistics of dimension nHypothesisxnHypothesis)
    if(test.type == "Dunnett" || test.type == "Parametric" || test.type == "Partly-Parametric"){
      for (epIDX in 1:nEps){
        if (EpType[[epIDX]] == "Continuous") {
          epSig <- Arms.std.dev[[epIDX]]
          sigma_0 <- epSig[1]
          sigma_trt <- epSig[-1]
          capLambda <- (sigma_0^2 + sigma_trt^2 / allocRatio[-1])^(-1)
        } else if (EpType[[epIDX]] == "Binary") {
          pi_c <- prop.ctr[[epIDX]]
          capLambda <- (1 / (pi_c * (1 - pi_c))) * (allocRatio[-1] / (1 + allocRatio[-1]))
        }

        InfoMatrix <- sapply(ctrSS, function(x) {
          x * capLambda
        }) # row=hypothesis, col=looks

        # Sigma/Corr for ith stage incremental data
        sigmaZ <- matrix(NA, nrow = length(capLambda), ncol = length(capLambda))

        for (l in 1:length(capLambda))
        {
          for (m in l:length(capLambda)) {
            if (EpType[[epIDX]] == "Continuous") {
              sigmaZ[l, m] <- sigmaZ[m, l] <- varCovZ(
                EpType = "Continuous", i1 = l, k1 = 1,
                i2 = m, k2 = 1,
                sigma_0 = sigma_0, sigma_trt = sigma_trt,
                ctrSS = ctrSS, trtSS = trtSS, InfoMatrix = InfoMatrix
              )
            } else if (EpType[[epIDX]] == "Binary") {
              sigmaZ[l, m] <- sigmaZ[m, l] <- varCovZ(
                EpType = "Binary", i1 = l, k1 = 1,
                i2 = m, k2 = 1,
                ctrProp = pi_c,
                ctrSS = ctrSS, trtSS = trtSS, InfoMatrix = InfoMatrix
              )
            }
          }
        }
        SigmaZ[[paste("EP", epIDX, sep = "")]] <- sigmaZ
      }
       Sigmalk <- as.matrix(Matrix::bdiag(SigmaZ))
       Sigmalk[Sigmalk == 0] <- NA
       rownames(Sigmalk) <- colnames(Sigmalk) <- paste("Z", 1:nrow(Sigmalk), sep = "")
      }else if(test.type == "Bonf"){
        Sigmalk <- matrix(NA, nrow = nHypothesis, ncol = nHypothesis)
        diag(Sigmalk) <- 1
        rownames(Sigmalk) <- colnames(Sigmalk) <- paste("Z", 1:nrow(Sigmalk), sep = "")
      }else if(test.type == "Sidak" || test.type == "Simes"){
        Sigmalk <- NA
      }
    Sigma[[paste("Stage", lkIDX, sep = "")]] <- Sigmalk
  }
  Sigma
}
