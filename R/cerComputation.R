# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

# Compute the CER for all intersection hypothesis after stage-1 analysis(only for CER analysis tool)
#' @param b2 stage-2 planned boundaries
#' @param WH Weights for all intersection hypothesis
#' @param p observed raw p-values of the current stage
#' @param test.type test type
#' @param HypoMap description for arms, endpoints and hypothesis
#' @param CommonStdDev Common Standard deviation flag
#' @param allocRatio planned allocation ratio
#' @param sigma planned arm-wise std.dev
#' @param Sigma planned var-cov
#' @param AllocSampleSize planned arm-wise sample size(cum.)
getCER <- function(b2,WH,p1,test.type,HypoMap,CommonStdDev,
                   allocRatio,sigma,Sigma,AllocSampleSize,EpType,prop.ctr){

  SUBSETS <- CONDERR <- c()
  Stage2Sigma <- getStage2PlanSigma(CommonStdDev = CommonStdDev,
                                    allocRatio = allocRatio,
                                    sigma = sigma,
                                    Sigma = Sigma,
                                    AllocSampleSize = AllocSampleSize,
                                    EpType = EpType,
                                    prop.ctr = prop.ctr)

  PlanSSHyp <- getHypoSS(SS = AllocSampleSize,
                         HypoMap = HypoMap)

  for(hypIDX in 1:nrow(WH)){
    n <- ncol(WH)
    J <- WH[hypIDX, 1:(n/2)]
    w1 <-  as.numeric(WH[hypIDX, (n/2+1):n])
    a2 <- as.numeric(b2[hypIDX,])

    get_Sets <- connSets(J = J, w = w1, test.type = test.type, HypoMap = HypoMap)
    conn_Sets <- get_Sets$connSets
    conn_Sets_name <- paste("(", paste(
      unlist(lapply(conn_Sets, function(x) {
        paste(x, collapse = ",")
      })),
      collapse = "),("
    ), ")", sep = "")

    SubSets <- conn_Sets_name

    testGrps <- connParamNParmSets(conn_Sets = conn_Sets)
    ParamGrps <- testGrps$ParamGrps
    NParamGrps <- testGrps$NParamGrp

    cerParamGrps <- pcerNParamGrps <- c()

    ################### CER for parametric subsets######################
    if (length(ParamGrps) != 0) {
      for (pGrp in ParamGrps) {
        if (length(pGrp) != 0) {
          Jh <- pGrp
          epIDX <- unique(HypoMap[pGrp, ]$Groups)
          wJh <- w1[pGrp]
          aJh <- a2[pGrp]
          cJ2Ratio <- aJh / wJh
          if (length(cJ2Ratio) > 1) {
            cJ2Ratio <- as.numeric(cJ2Ratio)
            if (max(abs(diff(cJ2Ratio))) > 1E-3) print("Error in Parametric SubGroup")
          }

          cJ2 <- cJ2Ratio[1]
          pJh <- p1[pGrp]

          # Compute Parametric CER based on old weights
          stage2sigmaS <- Stage2Sigma$SigmaSIncr[[epIDX]][floor(pGrp / epIDX), floor(pGrp / epIDX)]
          InfoMatrix <- Sigma$InfoMatrix[[epIDX]][floor(pGrp / epIDX), ]

          cerParam <- exitProbStage2Cond(
            cJ2 = cJ2, p1 = pJh, w = wJh,
            InfoMatrix = InfoMatrix, stage2sigmaS = stage2sigmaS, Conditional = TRUE
          )
          cerParamGrps <- c(cerParam, cerParamGrps)
        }
      }
    }
    ############# End of Parametric Computations ##################

    ############ PCER for non-parametric subsets######################

    if (length(NParamGrps) != 0) {
      wJh <- as.numeric(w1)[NParamGrps]
      aJh <- as.numeric(a2[NParamGrps])
      pJh <- as.numeric(p1[NParamGrps])

      pcer <- unlist(lapply(1:length(NParamGrps), function(x) {
        getPCER(a2 = aJh[x], p1 = pJh[x], ss1 = PlanSSHyp[[1]][x], ss2 = PlanSSHyp[[2]][x])
      }))

      pcerNParamGrps <- pcer
    }
    ############# End of Non-Parametric Computations ##################

    SubSets <- paste(
      paste("P :", paste(
        unlist(lapply(ParamGrps, function(x) {
          paste(x, collapse = ",")
        })),
        collapse = ","
      ), sep = ""),
      ",",
      paste("NP:", paste(
        unlist(lapply(NParamGrps, function(x) {
          paste(x, collapse = ",")
        })),
        collapse = ","
      ), sep = ""),
      sep = ""
    )
    roundDigit <- function(err, digits) {
      if (length(err) != 0) {
        round(err, digits)
      } else {
        0
      }
    }

    ConditionalError <- paste(
      paste("CER :", paste(
        unlist(lapply(roundDigit(cerParamGrps, 5), function(x) {
          paste(x, collapse = ",")
        })),
        collapse = ","
      ), sep = ""),
      ",",
      paste("PCER :", paste(
        unlist(lapply(roundDigit(pcerNParamGrps, 5), function(x) {
          paste(x, collapse = ",")
        })),
        collapse = ","
      ), sep = ""),
      sep = ""
    )

    SUBSETS <- c(SUBSETS, SubSets)
    CONDERR <- c(CONDERR, ConditionalError)
  }

  # Modified weights table
  HypoTab <- WH[, grep("H", names(WH))]

  InterHyp <- apply(HypoTab, 1, function(h) {
    paste(names(HypoTab)[which(h == 1)], collapse = ",")
  })

  InterWeight <- apply(WH, 1, function(h) {
    J <- which(h[1:(length(h) / 2)] == 1)
    w <- h[((length(h) / 2) + 1):length(h)]
    w <- sapply(w,function(x) roundDigit(x,3))
    paste(w[J], collapse = ",")
  })

  knitr::kable(data.frame("Hypotheses" = InterHyp,
             "Weights"=InterWeight,
             "SubSets"=SUBSETS,
             "CER"=CONDERR),
             align = "c")
}



getStage2PlanSigma <- function(CommonStdDev,allocRatio, sigma,
                               Sigma, AllocSampleSize, EpType, prop.ctr)
{
  #CommonStdDev flag is global variable set from simMAMSMEP(.), adaptGMCP_CER(.) function call
  if(CommonStdDev == T){
    for(sigIDX in 1:length(sigma))
      if(all(!is.na(sigma[[sigIDX]]))){
        sigma[[sigIDX]] <- sapply(sigma[[sigIDX]], function(x)sigma[[sigIDX]][1])
      }
  }
  ############################################

  SigmaZIncr <- list() # Z-Scale Incremental sigma
  SigmaSIncr <- list() # Score Incremental Scale sigma
  ############# With planned Samples for CER computations##############
  nGrps <- length(Sigma$SigmaZ)
  for (i in 1:nGrps) # For all groups
  {

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
  list("SigmaZIncr"=SigmaZIncr,
       "SigmaSIncr"=SigmaSIncr,
       "InfoMatrix"=InfoMatrix)
}



