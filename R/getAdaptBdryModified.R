# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

#Compute adaptive boundaries modified 12-06-2024
getAdaptBdryModified <- function(J, w1, w2, a2, a1, p1, test.type, HypoMap,
                         Sigma, Stage2Sigma, Stage2HypoIDX, PlanSSHyp, ModSSHyp) {

  ##Compute the CER and PCER based on Stage-1 design ##
  #get the stage-1 subsets from w1
  get_Sets_Stage1 <- connSets(J = J, w = w1, test.type = test.type, HypoMap = HypoMap)
  conn_Sets_Stage1 <- get_Sets_Stage1$connSets
  conn_Sets_Stage1_name <- paste("(", paste(
    unlist(lapply(conn_Sets_Stage1, function(x) {
      paste(x, collapse = ",")
    })),
    collapse = "),("
  ), ")", sep = "")
  SubSets_Stage1 <- conn_Sets_Stage1_name

  #Compute conditional errors CER/PCER
  CondError <- lapply(conn_Sets_Stage1, function(subSet){
    getCondErrorSubSets(subSet= subSet,HypoMap = HypoMap, w1 = w1,
                         a2 = a2,p1 = p1,Stage2Sigma=Stage2Sigma,PlanSSHyp=PlanSSHyp)
    })
  #----------------------------------------------------

  ##Compute the adaptive boundary based on Stage-2 design ##
  #get the stage-1 subsets from w1
  J_Plus <- rep(0,length(J))
  J_Plus[Stage2HypoIDX] <- 1
  get_Sets_Stage2 <- connSets(J = J_Plus, w = w2, test.type = test.type, HypoMap = HypoMap)
  conn_Sets_Stage2 <- get_Sets_Stage2$connSets
  conn_Sets_Stage2_name <- paste("(", paste(
    unlist(lapply(conn_Sets_Stage2, function(x) {
      paste(x, collapse = ",")
    })),
    collapse = "),("
  ), ")", sep = "")
  SubSets_Stage2 <- conn_Sets_Stage2_name
  Method <- get_Sets_Stage2$Method

  #two groups :  parametric(any conn_sets with length>1), non-parametric(any conn_sets with length = 1)
  testGrps <- connParamNParmSets(conn_Sets = conn_Sets_Stage2)
  ParamGrps <- testGrps$ParamGrps
  CondErrParam <- lapply(ParamGrps, function(pGrp){
    psetIDX <- vecMatch(pGrp, conn_Sets_Stage1)
    if(any(psetIDX)){
     return(CondError[[which(psetIDX)]])
    }else{
      return(NA)
    }
  })
  NParamGrps <- testGrps$NParamGrp
  CondErrNParam <- lapply(NParamGrps, function(npGrp){
    npsetIDX <- vecMatch(npGrp, conn_Sets_Stage1)
    if(any(npsetIDX)){
      return(CondError[[which(npsetIDX)]])
    }else{
      return(0)
    }
  })

  #Compute adjusted boundary

  Stage2AdjBdry <- rep(0, length(J))
  cerParamGrps <- pcerNParamGrps <- c()
  ScaleWeights <- c()

  ################### Adjusted Boundary for parametric subsets######################
  if (length(ParamGrps) != 0) {
    for (pGrp in ParamGrps) {
      if (length(pGrp) != 0) {
        Jh <- pGrp
        epIDX <- unique(HypoMap[pGrp, ]$Groups)
        wJh <- w1[pGrp]
        aJh <- a2[pGrp]
        cJ2Ratio <- aJh / wJh
        if (length(cJ2Ratio) > 1) {
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

        # Compute Stage-2 adaptive boundary based on new weights & distribution
        # boundaries for only available hypothesis which have non-zero weights
        pGrpMod <- pGrp[pGrp %in% Stage2HypoIDX &
                          pGrp %in% which(w2 > 0)]

        if (length(pGrpMod) > 1) {
          stage2sigmaSMod <- Stage2Sigma$Stage2SigmaS[[epIDX]][floor(pGrpMod / epIDX), floor(pGrpMod / epIDX)]

          InfoMatrixMod <- cbind(
            Sigma$InfoMatrix[[epIDX]][, 1],
            Stage2Sigma$Stage2InfoMatrixCum[[epIDX]]
          )[floor(pGrpMod / epIDX), ]

          cJ2Mod <- getStage2CondParamBdry(
            cer = cerParam, p1 = p1[pGrpMod],
            w = w2[pGrpMod], InfoMatrix = InfoMatrixMod,
            stage2sigmaS = stage2sigmaSMod, Conditional = TRUE
          )
          Stage2AdjBdry[pGrpMod] <- cJ2Mod * w2[pGrpMod]
        } else if (length(pGrpMod) == 1) # If after adaptation pGrpMod becomes a singleton set
        {
          SS1Mod <- ModSSHyp[[1]][pGrpMod]
          SS2Mod <- ModSSHyp[[2]][pGrpMod]

          SingleParmOut <- getStage2CondNParamBdry(
            a1 = a1[pGrpMod], p1 = p1[pGrpMod],
            v = w2[pGrpMod], BJ = cerParam, SS1 = SS1Mod, SS2 = SS2Mod
          )

          # Alternative approch to compute the boundary
          # alternative <- getStage2CondNParamBdry1Hypo(p1=p1[pGrpMod],v=w2[pGrpMod],BJ=cerParam,
          #                                             SS1=SS1Mod, SS2=SS2Mod)

          Stage2AdjBdry[pGrpMod] <- SingleParmOut$Stage2AdjBdry
          ScaleWeights <- c(ScaleWeights, SingleParmOut$adjWeights)
        }
      }
    }
  }
  ############# End of Parametric Computations ##################

  ############ Adjusted Boundary for non-parametric subsets######################

  if (length(NParamGrps) != 0) {
    wJh <- as.numeric(w1)[NParamGrps]
    aJh <- as.numeric(a2[NParamGrps])
    pJh <- as.numeric(p1[NParamGrps])

    pcer <- unlist(lapply(1:length(NParamGrps), function(x) {
      getPCER(a2 = aJh[x], p1 = pJh[x], ss1 = PlanSSHyp[[1]][x], ss2 = PlanSSHyp[[2]][x])
    }))

    pcerNParamGrps <- pcer
    cerNParam <- sum(pcer)

    # Compute Stage-2 adaptive boundary based on new Sample Size & weights

    NPGrpsMod <- NParamGrps[NParamGrps %in% Stage2HypoIDX &
                              NParamGrps %in% which(w2 > 0)] # boundaries for only available hypothesis

    if (length(NPGrpsMod) != 0) # & cerNParam < 1
    {
      SS1Mod <- ModSSHyp[[1]][NPGrpsMod]
      SS2Mod <- ModSSHyp[[2]][NPGrpsMod]

      nParmOut <- getStage2CondNParamBdry(
        a1 = a1[NPGrpsMod],
        p1 = p1[NPGrpsMod],
        v = w2[NPGrpsMod],
        BJ = cerNParam,
        SS1 = SS1Mod,
        SS2 = SS2Mod
      )
      # Alternative approch to compute the boundary
      # if(length(NPGrpsMod)==1){
      #   alternative <- getStage2CondNParamBdry1Hypo(p1=p1[NPGrpsMod],
      #                                               v=w2[NPGrpsMod],
      #                                               BJ=cerNParam,
      #                                               SS1=SS1Mod,
      #                                               SS2=SS2Mod)
      # }

      Stage2AdjBdry[NPGrpsMod] <- nParmOut$Stage2AdjBdry
      #ScaleWeights <- c(ScaleWeights, nParmOut$adjWeights)
    } else {
      Stage2AdjBdry[NPGrpsMod] <- 0
    }
  }


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
  if (length(ScaleWeights) != 0) {
    ScaleWeights2 <- paste0("(", paste0(roundDigit(ScaleWeights, 5), collapse = ","), ")")
  } else {
    ScaleWeights2 <- NA
  }


  list(
    "SubSets" = SubSets,
    "ConditionalError" = ConditionalError,
    "Stage2AdjBdry" = Stage2AdjBdry
    #"ScaledWeights" = ScaleWeights2
  )
}

getAdjBoundariesSubSets <- function(subSet_Stage2,conn_Sets_Stage1,CondError){

  CondError

}

vecMatch <- function(subSet,SET){
  unlist(lapply(SET, function(set){
    any(subSet %in% set)
  }))
}

getCondErrorSubSets <- function(subSet, HypoMap, w1, a2, p1, Stage2Sigma, PlanSSHyp){
  if(length(subSet) > 1){
    #Parametric CER
    pGrp <- subSet
    epIDX <- unique(HypoMap[pGrp, ]$Groups)
    wJh <- w1[pGrp]
    aJh <- a2[pGrp]
    cJ2Ratio <- aJh / wJh
    if (length(cJ2Ratio) > 1) {
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
    return(cerParam)
  }else{
    #Non-Parametric PCER
    NParamGrp <- subSet
    wJh <- as.numeric(w1)[NParamGrp]
    aJh <- as.numeric(a2[NParamGrp])
    pJh <- as.numeric(p1[NParamGrp])
    pcerNParam <- getPCER(a2 = aJh,p1 = pJh,
                    ss1 = PlanSSHyp[[1]][NParamGrp],
                    ss2 = PlanSSHyp[[2]][NParamGrp])
    return(pcerNParam)
  }
}


