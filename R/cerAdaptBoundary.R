#Boundary computation based on the mcpObj
adaptBdryCER <- function(mcpObj) {
  nHypothesis <- length(mcpObj$IntialHypothesis)
  nLooks <- length(mcpObj$Stage1Obj$info_frac)

  ## The weights before updating the testing strategy
  WH_old <- mcpObj$WH_Prev

  ## The weights after modifying the testing strategy with the targeted intersect hypotheses
  WH_modified <- mcpObj$WH
  Stage2PlanBdry <- mcpObj$Stage1Obj$plan_Bdry$Stage2Bdry
  Stage1PlanBdry <- mcpObj$Stage1Obj$plan_Bdry$Stage1Bdry

  WH_modified_idx <- as.vector(apply(WH_modified[, grep("H", names(WH_modified))], 1, function(x) {
    paste(x, collapse = "")
  }))
  WH_old_idx <- as.vector(apply(WH_old[, grep("H", names(WH_old))], 1, function(x) {
    paste(x, collapse = "")
  }))

  # Computed covariance matrix
  if (mcpObj$test.type == "Partly-Parametric" || mcpObj$test.type == "Parametric") {
    Stage2Sigma <- getStage2Sigma(
      nHypothesis = nHypothesis,
      EpType = mcpObj$lEpType,
      nLooks = nLooks,
      Sigma = mcpObj$Stage1Obj$Sigma,
      AllocSampleSize = mcpObj$AllocSampleSize,
      allocRatio = mcpObj$allocRatio,
      sigma = mcpObj$sigma,
      prop.ctr = mcpObj$prop.ctr,
      Stage2AllocSampleSize = mcpObj$Stage2AllocSampleSize,
      Stage2allocRatio = mcpObj$Stage2allocRatio,
      Stage2sigma = mcpObj$sigma,
      CommonStdDev = mcpObj$CommonStdDev
    )
  } else {
    Stage2Sigma <- NA
  }

  # Compute PCER based on planned sample size
  PlanSSHyp <- getHypoSS(SS = mcpObj$AllocSampleSize, HypoMap = mcpObj$HypoMap)
  ModSSHyp <- getHypoSS(SS = mcpObj$Stage2AllocSampleSize, HypoMap = mcpObj$HypoMap)

  SubSets <- c()
  ConditionalError <- c()
  Stage2AdjBdry <- c()
  ScaleWeights <- rep(NA, nrow(WH_modified))

  for (i in seq_len(nrow(WH_modified))) {
    # print(i)
    p1 <- as.numeric(mcpObj$p_raw)

    # Intersection hypothesis to test
    J <- as.numeric(WH_modified[i, grep("H", names(WH_modified))])
    oldIdx <- which(WH_old_idx == WH_modified_idx[i]) # Index for the WH_old corresponding to J
    w1 <- as.numeric(WH_old[oldIdx, grep("Weight", names(WH_old))])
    a2 <- as.numeric(Stage2PlanBdry[oldIdx, ])
    a1 <- as.numeric(Stage1PlanBdry[oldIdx, ])

    # Available hypothesis for stage-2
    I2 <- get_numeric_part(mcpObj$IndexSet)
    J2 <- rep(0, length(J))
    J2[I2] <- 1
    J2 <- J * J2 # Weights for the intersection J,J2

    modIdx <- which(WH_modified_idx == paste(J2, collapse = "")) # Index for the WH_old corresponding to J
    w2 <- as.numeric(WH_modified[modIdx, grep("Weight", names(WH_modified))])



    adaptOut <- getAdaptBdry(
      J = J, w1 = w1, w2 = w2, a2 = a2, a1 = a1, p1 = p1,
      test.type = mcpObj$test.type,
      HypoMap = mcpObj$HypoMap,
      Sigma = mcpObj$Stage1Obj$Sigma,
      Stage2Sigma = Stage2Sigma,
      PlanSSHyp = PlanSSHyp,
      ModSSHyp = ModSSHyp,
      Stage2HypoIDX = get_numeric_part(mcpObj$IndexSet)
    )

    SubSets <- c(SubSets, adaptOut$SubSets)
    ConditionalError <- c(ConditionalError, adaptOut$ConditionalError)
    Stage2AdjBdry <- rbind(Stage2AdjBdry, adaptOut$Stage2AdjBdry)
    ScaleWeights[i] <- adaptOut$ScaledWeights
  }

  colnames(Stage2AdjBdry) <- paste("a", 1:nHypothesis, "2_adj", sep = "")

  # Modified Sample Size table

  ModfiedSSTab <- knitr::kable(mcpObj$Stage2AllocSampleSize, align = "c")

  # Modified weights table
  HypoTab <- WH_modified[, grep("H", names(WH_modified))]

  InterHyp <- apply(HypoTab, 1, function(h) {
    paste(names(HypoTab)[which(h == 1)], collapse = ",")
  })

  InterWeight <- apply(WH_modified, 1, function(h) {
    J <- which(h[1:(length(h) / 2)] == 1)
    w <- h[((length(h) / 2) + 1):length(h)]
    paste(w[J], collapse = ",")
  })

  ModifiedWeightTab <- WH_modified


  # Adaptive Test Procedure
  # TestProcedureTab <- cbind(WH_modified[1:nHypothesis],SubSets,ConditionalError)
  TestProcedureTab1 <- knitr::kable(data.frame(
    "Hypotheses" = InterHyp,
    "Weights" = InterWeight,
    "ScaleWeights" = ScaleWeights,
    "SubSets" = SubSets,
    "Conditional_Error" = ConditionalError,
    row.names = NULL
  ), align = "c")


  # Adjusted Boundary
  # AdjBdryTab <- cbind(WH_modified[1:nHypothesis],Stage2AdjBdry)
  stg2bdry <- c()

  for (hypIDX in 1:nrow(Stage2AdjBdry)) {
    J <- which(WH_modified[hypIDX, 1:(length(WH_modified) / 2)] == 1)
    stg2bdry <- c(
      stg2bdry,
      paste(round(Stage2AdjBdry[hypIDX, J], 6), collapse = ",")
    )
  }
  AdjBdryTab1 <- knitr::kable(data.frame(
    "Hypotheses" = InterHyp,
    "Adapt_Boundary" = stg2bdry,
    row.names = NULL
  ), align = "c")


  AdaptTable <- list(
    "Sample_Size" = ModfiedSSTab,
    "Stage2_Test_Procedure" = TestProcedureTab1,
    "Adjusted_Boundary" = AdjBdryTab1
  )

  list(
    "Stage2Tables" = AdaptTable,
    "Stage2AdjBdry" = Stage2AdjBdry,
    "Stage2Sigma" = Stage2Sigma
  )
}




getAdaptBdry <- function(J, w1, w2, a2, a1, p1, test.type, HypoMap,
                         Sigma, Stage2Sigma, Stage2HypoIDX, PlanSSHyp, ModSSHyp) {
  get_Sets <- connSets(J = J, w = w1, test.type = test.type, HypoMap = HypoMap)
  conn_Sets <- get_Sets$connSets
  conn_Sets_name <- paste("(", paste(
    unlist(lapply(conn_Sets, function(x) {
      paste(x, collapse = ",")
    })),
    collapse = "),("
  ), ")", sep = "")

  SubSets <- conn_Sets_name
  Method <- get_Sets$Method

  testGrps <- connParamNParmSets(conn_Sets = conn_Sets)
  ParamGrps <- testGrps$ParamGrps
  NParamGrps <- testGrps$NParamGrp

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
        ScaleWeights <- c(ScaleWeights, nParmOut$adjWeights)
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
    "Stage2AdjBdry" = Stage2AdjBdry,
    "ScaledWeights" = ScaleWeights2
  )
}


# to get samples hypothesis wise
getHypoSS <- function(SS, HypoMap) {
  SS_H <- lapply(1:nrow(SS), function(j) {
    unlist(lapply(1:nrow(HypoMap), function(i) {
      sum(SS[j, ][as.numeric(HypoMap[i, c("Control", "Treatment")])], na.rm = T)
    }))
  })
  SS_H
}
