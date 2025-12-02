# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

# The file contains supporting functions for Adaptive GMCP Simulation ----
## Author: Ajoy.M
getPreSimObjs <- function(gmcpSimObj) {
  #----------------------------------------------------------------------------
  ################## Generate unique Seed for the run #########################
  #should not be greater than 4 digits
  if(gmcpSimObj$Seed == "Random"){
   IntSeed1 <- as.integer(Sys.time())
   IntSeed2 <- ifelse(IntSeed1 > 9999,
        substr(as.character(IntSeed1), nchar(IntSeed1) - 3, nchar(IntSeed1)),
        substr(IntSeed1)
        )
   SimSeed <-as.integer(
                paste(sample(strsplit(IntSeed2, "")[[1]]),collapse = '')
                )

  } else {
    SimSeed <- gmcpSimObj$Seed
  }

  #----------------------------------------------------------------------------
  ################### Arms-EPs-Hypotheses Mapping #####################################
  HypoMap <- getHypoMap2(
    des.type = gmcpSimObj$des.type,
    nHypothesis = gmcpSimObj$nHypothesis,
    nEps = gmcpSimObj$nEps,
    nArms = gmcpSimObj$nArms,
    lEpType = gmcpSimObj$lEpType
  )


  #------------------------------------------------------------------------------
  #################### Identify True Null based on given Response #######################
  # TrueNull <- checkTrueNull2(HypoMap = HypoMap, Arms.Mean = gmcpSimObj$Arms.Mean)
  # returns a vector of boolean of length nHypothesis indicating whether the hypothesis is true null
  TrueNull <- checkTrueNull3(
    HypoMap = HypoMap,
    Arms.Mean = gmcpSimObj$Arms.Mean,
    Arms.Prop = gmcpSimObj$Arms.Prop
  )


  #---------------------------------------------------------------------------------------
  #################### Weights for all intersection hypothesis ##########################
  ############# Weights for all intersection hypothesis#########################
  allGraphs <- genWeights(
    w = gmcpSimObj$IntialWeights,
    g = gmcpSimObj$G,
    HypothesisName = HypoMap$Hypothesis
  )
  # This is a dataframe containing all the intersection hypothesis and their weights
  WH <- allGraphs$IntersectionWeights

  # WH <- as.data.frame(gMCPLite::generateWeights(g = gmcpSimObj$G,w = gmcpSimObj$IntialWeights))
  # colnames(WH) <- c(paste("H",1:gmcpSimObj$nHypothesis, sep=''),
  #                  paste("Weight",1:gmcpSimObj$nHypothesis, sep = ''))

  #--------------------------------------------------------------------------------------
  ############################## Planned Sample Size ################################

  planSS <- getPlanAllocatedSamples(
    SS = gmcpSimObj$Max_SS,
    allocRatio = gmcpSimObj$Arms.alloc.ratio,
    info_frac = gmcpSimObj$InfoFrac
  )

  planSSIncr <- planSS$IncrementalSamples
  planSSCum <- planSS$CumulativeSamples

  #----------------------------------------------------------------------------------
  ######## Computation of Inverse Normal Weights from planned sample size##############
  # Needed for combining p-values method and CER with FWER control as 'CombinationTest'
  InvNormWeights <- getInvNormWeights(planSSIncr = planSSIncr)

  #----------------------------------------------------------------------------------

  # Prepare the preSim objects for two methods
  if (gmcpSimObj$Method == "CombPValue") # Inputs for Combining P-values method
    {
      #----------------------------------------------------------------------------------
      ###################### Get the stage-wise p-value boundaries####################
      pValBdry <- getPvalBdry(
        alpha = gmcpSimObj$alpha,
        nLooks = gmcpSimObj$nLooks,
        info_frac = gmcpSimObj$InfoFrac,
        typeOfDesign = gmcpSimObj$typeOfDesign,
        deltaWT = gmcpSimObj$deltaWT,
        deltaPT1 = gmcpSimObj$deltaPT1,
        gammaA = gmcpSimObj$gammaA,
        userAlphaSpending = gmcpSimObj$userAlphaSpending
      )

      #----------------------------------------------------------------------------------
      ###################### Get the stage-wise incremental Correlation####################
      PlanCorrelation <- getPlanCorrelation(
        nHypothesis = gmcpSimObj$nHypothesis,
        SS_Incr = planSSIncr,
        Arms.std.dev = gmcpSimObj$Arms.std.dev,
        test.type = gmcpSimObj$test.type,
        EpType = gmcpSimObj$lEpType,
        prop.ctr = gmcpSimObj$prop.ctr,
        CommonStdDev = gmcpSimObj$CommonStdDev
      )

      #----------------------------------------------------------------------------------
      PreSimObj <- list(
        "IndexSet" = HypoMap$Hypothesis,
        "SimSeed" = SimSeed, "HypoMap" = HypoMap,
        "TrueNull" = TrueNull, "planSS" = planSS,
        "WH" = WH, "InvNormWeights" = InvNormWeights,
        "pValBdry" = pValBdry, "PlanCorrelation" = PlanCorrelation,
        "W_Norm" = InvNormWeights$W_Norm,
        "AdjPValues" = NA
      )
    } else if (gmcpSimObj$Method == "CER") # Inputs for CER method
    {
      #----------------------------------------------------------------------------------
      ######## Computation of covariance matrix##############
      # Bug fix: Earlier Sigma was calculated only in case of Partly-Parametric
      # and Parametric tests and it was set to NA in case of Non-Parametric.
      # However, this led to problems in the function SingleSimCER() as it
      # requires the fisher info matrix in all cases.
      # Fixed this bug by calculating Sigma unconditionally.
      # if (gmcpSimObj$test.type == "Partly-Parametric" || gmcpSimObj$test.type == "Parametric") {
        Sigma <- getSigma(
          SS_Cum = planSS$CumulativeSamples,
          sigma = gmcpSimObj$Arms.std.dev,
          allocRatio = gmcpSimObj$Arms.alloc.ratio,
          EpType = gmcpSimObj$lEpType,
          prop.ctr = gmcpSimObj$prop.ctr,
          CommonStdDev = gmcpSimObj$CommonStdDev,
          info_frac = gmcpSimObj$InfoFrac
        )
      # } else {
      #   Sigma <- NA
      # }
      #----------------------------------------------------------------------------------
      ######## Computation of Planned Boundaries##############
      plan_Bdry <- planBdryCER(
        nHypothesis = gmcpSimObj$nHypothesis,
        nEps = gmcpSimObj$nEps,
        nLooks = gmcpSimObj$nLooks,
        alpha = gmcpSimObj$alpha,
        info_frac = gmcpSimObj$InfoFrac,
        typeOfDesign = gmcpSimObj$typeOfDesign,
        deltaWT = gmcpSimObj$deltaWT,
        deltaPT1 = gmcpSimObj$deltaPT1,
        gammaA = gmcpSimObj$gammaA,
        userAlphaSpending = gmcpSimObj$userAlphaSpending,
        test.type = gmcpSimObj$test.type,
        Sigma = Sigma,
        WH = WH,
        HypoMap = HypoMap,
        Scale = "Score",
        planSSCum = planSSCum
      )

      # ###################################
      # # Ani: hardcoding the following variables for debugging
      # plan_Bdry$Stage2Bdry[1,1] <- plan_Bdry$Stage2Bdry[1,2] <- 0.0131202663476177
      # plan_Bdry$Stage2Bdry[2,2] <- plan_Bdry$Stage2Bdry[3,1] <- 0.0244997717772641
      # ###################################

      #----------------------------------------------------------------------------------
      PreSimObj <- list(
        "IndexSet" = HypoMap$Hypothesis,
        "IntialHypothesis" = HypoMap$Hypothesis,
        "SimSeed" = SimSeed, "HypoMap" = HypoMap,
        "TrueNull" = TrueNull, "planSS" = planSS,
        "WH" = WH, "Sigma" = Sigma,
        "plan_Bdry" = plan_Bdry, "InvNormWeights" = InvNormWeights
      )
    }

  # Plot the initial Graph
  SubText <- getPlotText(HypoMap)
  if (gmcpSimObj$plotGraphs) {
    iniGraph <- plotGraph(
      HypothesisName = HypoMap$Hypothesis,
      w = gmcpSimObj$IntialWeights,
      G = gmcpSimObj$G,
      Title = "Initial Graph",
      Text = SubText
    )
    PreSimObj$iniGraph <- iniGraph
  }

  PreSimObj
}
