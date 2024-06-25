#### Compute plan boundary for parametric, non-parametric and mixed case
planBdryCER <- function(nHypothesis,
                        nEps,
                        nLooks,
                        alpha,
                        info_frac,
                        typeOfDesign,
                        deltaWT,
                        deltaPT1,
                        gammaA,
                        userAlphaSpending,
                        test.type,
                        Sigma,
                        WH,
                        HypoMap,
                        Scale) {
  Method <- c()
  SubSets <- c()
  SignfLevel <- c()

  Stage1Bdry <- Stage2Bdry <- matrix(0, nrow = nrow(WH), ncol = nHypothesis)
  colnames(Stage1Bdry) <- paste("a", 1:nHypothesis, "1", sep = "")
  colnames(Stage2Bdry) <- paste("a", 1:nHypothesis, "2", sep = "")

  for (i in 1:nrow(WH))
  {
    # print(i)
    J <- as.numeric(WH[i, 1:nHypothesis])
    w <- as.numeric(WH[i, (nHypothesis + 1):(2 * nHypothesis)])

    get_Sets <- connSets(J = J, w = w, test.type = test.type, HypoMap = HypoMap)
    conn_Sets <- get_Sets$connSets
    conn_Sets_name <- paste("(", paste(
      unlist(lapply(conn_Sets, function(x) {
        paste(x, collapse = ",")
      })),
      collapse = "),("
    ), ")", sep = "")

    SubSets <- c(SubSets, conn_Sets_name)
    Method <- c(Method, get_Sets$Method)

    siglev <- list()
    for (edx in conn_Sets) {
      if (length(edx) > 1) {
        ## Dunnett Weighted Parametric ##
        sig_level <- alpha * sum(w[edx])
        gIDX <- unique(HypoMap[edx, ]$Groups)
        wJ <- as.numeric(w) # wJ : weights for the intersection J(including weights with 0)

        plan_parm_bdry <- getPlanParmBdry(
          gIDX = gIDX, hIDX = edx, alpha = sig_level, nLooks = nLooks,
          info_frac = info_frac, wJ = wJ,
          Sigma = Sigma, typeOfDesign = typeOfDesign,deltaWT= deltaWT,
          deltaPT1 = deltaPT1, gammaA = gammaA,
          userAlphaSpending = userAlphaSpending, Scale = Scale
        )
        Stage1Bdry[i, edx] <- plan_parm_bdry$Stage1Bdry[edx]
        Stage2Bdry[i, edx] <- plan_parm_bdry$Stage2Bdry[edx]
      } else {
        ## Non-Parametric ##
        sig_level <- alpha * sum(w[edx])
        plan_nparm_bdry <- getPlanNonParmBdry(
          nLooks = nLooks, sig_level = sig_level,
          info_frac = info_frac, typeOfDesign = typeOfDesign, deltaWT= deltaWT,
          deltaPT1 = deltaPT1, gammaA = gammaA,
          userAlphaSpending = userAlphaSpending
        )
        Stage1Bdry[i, edx] <- plan_nparm_bdry$Stage1Bdry
        Stage2Bdry[i, edx] <- plan_nparm_bdry$Stage2Bdry
      }
      siglev <- append(siglev, sig_level)
    }
    SignfLevel <- c(
      SignfLevel,
      paste("(", paste(
        unlist(lapply(siglev, function(x) {
          paste(x, collapse = ",")
        })),
        collapse = "),("
      ), ")", sep = "")
    )
  }
  #### Preparation of output tables ###

  # Intersection Weights#
  InterWeightsTab <- WH

  # Weight Table Preparation for outputs
  HypoTab <- WH[, 1:(ncol(WH) / 2)]
  InterHyp <- apply(HypoTab, 1, function(h) {
    paste(names(HypoTab)[which(h == 1)], collapse = ",")
  })
  InterWeight <- apply(WH, 1, function(h) {
    J <- which(h[1:(length(h) / 2)] == 1)
    w <- h[((length(h) / 2) + 1):length(h)]
    paste(w[J], collapse = ",")
  })
  WeightTab <- data.frame(
    "Hypothesis" = InterHyp,
    "Weights" = InterWeight,
    row.names = NULL
  )


  # Test procedure
  TestProcedureTab <- cbind(
    WH[, 1:(ncol(WH) / 2)],
    SubSets,
    Method,
    SignfLevel
  )

  TestProcedureTab1 <- data.frame(
    "Hypotheses" = InterHyp, "Weights" = InterWeight,
    "SubSets" = SubSets, "Method" = Method,
    row.names = NULL
  )
  TestProcedureTab1 <- knitr::kable(TestProcedureTab1, align = "c")
  # Stage-1 Boundary
  Stage1BdryTab <- cbind(WH[, 1:(ncol(WH) / 2)], Stage1Bdry)
  stg1bdry <- c()

  for (hypIDX in 1:nrow(Stage1Bdry)) {
    J <- which(WH[hypIDX, 1:(length(WH) / 2)] == 1)
    stg1bdry <- c(
      stg1bdry,
      paste(round(Stage1Bdry[hypIDX, J], 6), collapse = ",")
    )
  }

  Stage1BdryTab1 <- knitr::kable(data.frame(
    "Hypotheses" = InterHyp,
    "Alpha" = SignfLevel,
    "Stage1Boundary" = stg1bdry,
    row.names = NULL
  ), align = "c")


  PlanBdryTab <- list(
    "Test_Procedure" = TestProcedureTab1,
    "Stage1_Bounday" = Stage1BdryTab1
  )
  if (nLooks == 2) {
    Stage2BdryTab <- cbind(WH[, 1:(ncol(WH) / 2)], Stage2Bdry)
    stg2bdry <- c()
    for (hypIDX in 1:nrow(Stage2Bdry)) {
      J <- which(WH[hypIDX, 1:(length(WH) / 2)] == 1)
      stg2bdry <- c(
        stg2bdry,
        paste(round(Stage2Bdry[hypIDX, J], 6), collapse = ",")
      )
    }
    Stage2BdryTab1 <- knitr::kable(data.frame(
      "Hypotheses" = InterHyp,
      "Alpha" = SignfLevel,
      "Stage2Boundary" = stg2bdry,
      row.names = NULL
    ), align = "c")

    PlanBdryTab <- list(
      "Test_Procedure" = TestProcedureTab1,
      "Stage1_Bounday" = Stage1BdryTab1,
      "Stage2_Bounday" = Stage2BdryTab1
    )
  }

  list(
    "Stage1Bdry" = Stage1Bdry,
    "Stage2Bdry" = Stage2Bdry,
    "PlanBdryTable" = PlanBdryTab
  )
}
