# The file contains supporting functions for Adaptive GMCP Simulation Function AdaptGmcp_Simulation/adaptGMCP_SIM(.) ----
## Author: Ajoy.M


#------------------ -
# Perform  iterative simulation for Multi-Arm designs
#------------------ -
mnMAMSMEP_sim <- function(gmcpSimObj) {
  SimSeed <- ifelse(gmcpSimObj$Seed == "Random", sample(1:10000, 1), gmcpSimObj$Seed)

  ############################## Planned Sample Size ################################
  plannedSS <- getPlanSS(gmcpSimObj)
  D <- length(gmcpSimObj$IntialWeights)
  K <- gmcpSimObj$nLooks

  ######## Computation of Inverse Normal Weights from planned sample size##############
  SS_looks_incr <- rowSums(plannedSS)
  SS_looks_cumulative <- cumsum(SS_looks_incr)
  W_Norm <- matrix(NA, nrow = K, ncol = K)

  for (i in 1:nrow(W_Norm))
  {
    for (j in 1:i) W_Norm[i, j] <- sqrt(SS_looks_incr[j] / SS_looks_cumulative[i])
  }
  InvNormWeights <- W_Norm
  colnames(InvNormWeights) <- paste("W", 1:K, sep = "")
  rownames(InvNormWeights) <- paste("Look", 1:K, sep = "")
  W_Norm <- W_Norm[-1, ] # Removing the first row

  ################################################################################

  ###################### Get the stage-wise p-value boundaries####################
  UseExternal <- T
  if (UseExternal) # this part of the code can be replaced later with the internal R-codes
    {
    if(typeOfDesign == "WT"){
      des <- rpact::getDesignGroupSequential(
        kMax = K, alpha = gmcpSimObj$alpha,
        informationRates = info_frac,
        typeOfDesign = typeOfDesign,
        deltaWT = deltaWT
      )
    }else if(typeOfDesign == "PT"){
      des <- rpact::getDesignGroupSequential(
        kMax = K, alpha = gmcpSimObj$alpha,
        informationRates = info_frac,
        typeOfDesign = typeOfDesign,
        deltaPT1 = deltaPT1
      )
    }else if(typeOfDesign == "asHSD" || typeOfDesign == "asKD"){
      des <- rpact::getDesignGroupSequential(
        kMax = K, alpha = gmcpSimObj$alpha,
        informationRates = info_frac,
        typeOfDesign = typeOfDesign,
        gammaA = gammaA
      )
    }else if(typeOfDesign == "asUser"){
      des <- rpact::getDesignGroupSequential(
        kMax = K, alpha = gmcpSimObj$alpha,
        informationRates = info_frac,
        typeOfDesign = typeOfDesign,
        userAlphaSpending = userAlphaSpending
      )
    }else{
      des <- rpact::getDesignGroupSequential(
        kMax = K, alpha = gmcpSimObj$alpha,
        informationRates = info_frac,
        typeOfDesign = typeOfDesign
      )
    }

    Threshold <- des$stageLevels
    incr_alpha <- c(des$alphaSpent[1], diff(des$alphaSpent))
    # BdryTab
    bdryTab <- data.frame(
      "Look" = 1:K, "Information_Fraction" = info_frac,
      "Incr_alpha_spent" = incr_alpha,
      "ZScale_Eff_Bbry" = des$criticalValues,
      "PValue_Eff_Bbry" = Threshold,
      row.names = NULL
    )
  }
  gmcpSimObj$EffCutOff <- Threshold

  ###############################################################################
  GlobalIndexSet <- paste("H", 1:D, sep = "")
  ## Hypothesis
  if (!is.list(gmcpSimObj$Hypothesis)) # Hypothesis == 'CommonControl'
    {
      pairs.test <- lapply(names(gmcpSimObj$Arms.Mean)[-1], function(x) {
        c(x, names(gmcpSimObj$Arms.Mean)[1])
      })
    } else {
    pairs.test <- gmcpSimObj$Hypothesis
  }
  index_map <- data.frame(
    "GlobalIndexSet" = GlobalIndexSet,
    "Names" = unlist(lapply(pairs.test, paste, collapse = "-")),
    row.names = NULL
  )



  #################### Weights for all intersection hypothesis ##########################
  ############# Weights for all intersection hypothesis#########################
  allGraphs <- genWeights(w = WI, g = G, HypothesisName = GlobalIndexSet)
  WH <- allGraphs$IntersectionWeights

  # WH <- as.data.frame(gMCPLite::generateWeights(g = gmcpSimObj$G,w = gmcpSimObj$IntialWeights))
  # colnames(WH) <- c(GlobalIndexSet, paste("Weight",1:D, sep = ''))

  TrueNullDf1 <- checkTrueNull(gmcpSimObj, index_map)
  TrueNull <- TrueNullDf1$TrueNull
  names(TrueNull) <- paste("RejStatus", get_numeric_part(TrueNullDf1$GlobalIndexSet), sep = "")


  PreSimObj <- list(
    "PlannedSS" = plannedSS, "WH" = WH, "W_Norm" = W_Norm,
    "GlobalIndexSet" = GlobalIndexSet, "index_map" = index_map,
    "TrueNull" = TrueNull
  )

  ################################ Simulation ######################################
  SummaryStatFile <- ArmWiseSummary <- PowerTab <- data.frame()
  powersName <- c("simID", "nG", "nC", "nD", "nF")
  PowerTab <- data.frame(matrix(nrow = 0, ncol = length(powersName)))
  names(PowerTab) <- powersName

  starttime <- Sys.time()
  out <- lapply(1:gmcpSimObj$nSimulation, function(x) {
    SingleSimAnalysis(x, gmcpSimObj, PreSimObj, SimSeed)
  })
  for (i in 1:length(out))
  {
    SummaryStatFile <- rbind(SummaryStatFile, out[[i]]$SummStatDF)
    ArmWiseSummary <- rbind(ArmWiseSummary, out[[i]]$ArmWiseDF)
    PowerTab <- rbind(PowerTab, out[[i]]$powerCountDF)
  }
  elapsedTime <- Sys.time() - starttime
  Sim_power <- SimPowers(nSimulation = gmcpSimObj$nSimulation, PowerTab = PowerTab)

  ########################### End of the Simulation ################################

  if (gmcpSimObj$SummaryStat) {
    return(
      list(
        "Boundary_Table" = bdryTab,
        "Inverse_Normal_Weights" = InvNormWeights,
        "Overall_Powers" = Sim_power,
        "Summary_Stat" = SummaryStatFile,
        "ArmWiseSummary" = ArmWiseSummary,
        "Seed" = SimSeed,
        "elapsedTime" = elapsedTime
      )
    )
  } else {
    return(
      list(
        "Boundary_Table" = bdryTab,
        "Inverse_Normal_Weights" = InvNormWeights,
        "Overall_Powers" = Sim_power,
        "Seed" = SimSeed,
        "elapsedTime" = elapsedTime
      )
    )
  }
}

#----------------- -
# Perform Single Simulation
#----------------- -
SingleSimAnalysis <- function(simID, gmcpSimObj, PreSimObj, SimSeed) {
  # print(simID)
  D <- length(gmcpSimObj$IntialWeights)
  K <- gmcpSimObj$nLooks
  ArmsIndex <- names(gmcpSimObj$Arms.Mean)

  ArmsPresent <- rep(T, gmcpSimObj$nArms)
  ArmsRetained <- rep(F, gmcpSimObj$nArms)
  names(ArmsPresent) <- names(ArmsRetained) <- ArmsIndex

  rej_flag_Prev <- rej_flag_Curr <- DropedFlag <- rep(FALSE, D)
  names(rej_flag_Prev) <- names(rej_flag_Curr) <- names(DropedFlag) <- paste("H", 1:D, sep = "")

  # Output data
  SummStatNames <- c(
    "SimID", "LookID",
    paste("TestStat", get_numeric_part(PreSimObj$GlobalIndexSet), sep = ""),
    paste("RawPvalues", get_numeric_part(PreSimObj$GlobalIndexSet), sep = ""),
    paste("RejStatus", get_numeric_part(PreSimObj$GlobalIndexSet), sep = "")
  )
  SummStatDF <- data.frame(matrix(nrow = 0, ncol = length(SummStatNames)))
  names(SummStatDF) <- SummStatNames

  ArmWiseDataName <- c("SimID", "LookID", "Arm", "Completers", "Mean", "SE")
  ArmWiseDF <- data.frame(matrix(nrow = 0, ncol = length(ArmWiseDataName)))
  names(ArmWiseDF) <- ArmWiseDataName

  # info to run per look analysis
  mcpObj <- list(
    "CurrentLook" = 1,
    "IntialWeights" = gmcpSimObj$IntialWeights,
    "IntialHypothesis" = PreSimObj$GlobalIndexSet,
    "test.type" = gmcpSimObj$test.type,
    "IndexSet" = PreSimObj$GlobalIndexSet,
    "ArmsPresent" = ArmsPresent,
    "ArmsRetained" = ArmsRetained,
    "ArmsIndex" = ArmsIndex,
    "p_raw" = NA,
    "WH_Prev" = PreSimObj$WH,
    "WH" = PreSimObj$WH,
    "Correlation" = gmcpSimObj$Correlation,
    "AdjPValues" = NA,
    "W_Norm" = PreSimObj$W_Norm,
    "CutOff" = NA,
    "MultipleWinners" = gmcpSimObj$MultipleWinners,
    "rej_flag_Prev" = rej_flag_Prev,
    "rej_flag_Curr" = rej_flag_Curr,
    "SelectionLook" = c(),
    "SelectedIndex" = NA,
    "DropedFlag" = DropedFlag,
    "LastLook" = gmcpSimObj$nLooks,
    "Modify" = F,
    "ModificationLook" = c(),
    "newWeights" = NA,
    "newG" = NA,
    "PlannedSS" = PreSimObj$PlannedSS,
    "SummStatNames" = SummStatNames,
    "SummStat" = NA,
    "ArmData" = NA,
    "ContTrial" = T
  )
  while (mcpObj$ContTrial) {
    # print(mcpObj$CurrentLook)
    LookData <- genLookData(simID = simID, mcpObj = mcpObj, PreSimObj = PreSimObj, gmcpSimObj = gmcpSimObj, SimSeed = SimSeed)
    SummStat <- LookData$SummStat
    ArmData <- LookData$ArmData
    mcpObj$SummStat <- SummStat
    mcpObj$ArmData <- ArmData

    p_raw <- unlist(SummStat[, grep("RawPvalues", names(SummStat))]) ## Simulated raw p-values
    names(p_raw) <- paste("H", get_numeric_part(names(p_raw)), sep = "")

    mcpObj$p_raw <- addNAPvalue(p_raw, PreSimObj$GlobalIndexSet)
    mcpObj$CutOff <- gmcpSimObj$EffCutOff[mcpObj$CurrentLook]

    mcpObj <- PerLookMCPAnalysis(mcpObj) # Perform Testing

    rejStatus <- data.frame(matrix(mcpObj$rej_flag_Curr, nrow = 1))
    names(rejStatus) <- paste("RejStatus", get_numeric_part(names(mcpObj$rej_flag_Curr)), sep = "")
    SummStat <- fillNa(1, SummStat, rejStatus)

    SummStatDF <- plyr::rbind.fill(SummStatDF, SummStat)
    ArmWiseDF <- plyr::rbind.fill(ArmWiseDF, ArmData)

    mcpObj$rej_flag_Prev <- mcpObj$rej_flag_Curr
    mcpObj$WH_Prev <- mcpObj$WH

    if (StopTrial(mcpObj)) {
      mcpObj$ContTrial <- F
      break # Early Stopping
    } else # Next look
    {
      # Available Arms after rejection
      mcpObj$ArmsPresent <- getDroppedArms(ArmsPresent = mcpObj$ArmsPresent, rejflags = mcpObj$rej_flag_Curr, index_map = PreSimObj$index_map)
      # Selection for next look
      if (gmcpSimObj$Selection & (mcpObj$CurrentLook < K)) {
        mcpObj <- do_SelectionSim(mcpObj, gmcpSimObj, PreSimObj)
      }
      # Modify the weights and testing strategy
      if (gmcpSimObj$UpdateStrategy & (mcpObj$CurrentLook < K) & (length(mcpObj$IndexSet) > 1)) {
        mcpObj <- do_modifyStrategy(mcpObj)
      }
    }
    mcpObj$CurrentLook <- mcpObj$CurrentLook + 1

    # End of Look-wise While Loop #
  }
  # Count Power(Efficacy)
  powerCountDF <- CountPower(simID = simID, SummaryStatFile = SummStatDF, TrueNull = PreSimObj$TrueNull)
  list("SummStatDF" = SummStatDF, "ArmWiseDF" = ArmWiseDF, "powerCountDF" = powerCountDF)
}

#------------------ -
# Replace row-specific NA values from one data frame(df1) with the values from another(df2)
#------------------ -
fillNa <- function(rowID, df1, df2) {
  commonCol <- intersect(names(df1), names(df2))
  df1[rowID, commonCol] <- df2[rowID, commonCol]
  df1
}

#------------------ -
# Generate Normal Response for each look
#------------------ -
genLookData <- function(simID, mcpObj, PreSimObj, gmcpSimObj, SimSeed) {
  lookID <- mcpObj$CurrentLook
  df <- perArmData(SimSeed, simID, lookID, mcpObj, gmcpSimObj)

  HypoToTest <- PreSimObj$index_map[PreSimObj$index_map$GlobalIndexSet %in% mcpObj$IndexSet, ]

  testStat <- unlist(lapply(HypoToTest$Names, function(x) {
    pairTestStat(x, df)
  }))
  testStat <- data.frame(matrix(testStat, nrow = 1))
  names(testStat) <- paste("TestStat", get_numeric_part(mcpObj$IndexSet), sep = "")

  is.Right <- gmcpSimObj$TailType == "RightTail"
  p_values <- unlist(lapply(testStat, pnorm, lower.tail = !is.Right))
  p_values <- data.frame(matrix(p_values, nrow = 1))
  names(p_values) <- paste("RawPvalues", get_numeric_part(mcpObj$IndexSet), sep = "")

  testStat <- cbind(
    data.frame("SimID" = simID, "LookID" = lookID),
    testStat, p_values
  )

  ArmData <- cbind(data.frame("SimID" = rep(simID, nrow(df)), "LookID" = rep(lookID, nrow(df))), df)

  SummStatNames <- mcpObj$SummStatNames
  SummStat <- data.frame(matrix(nrow = 0, ncol = length(SummStatNames)))
  names(SummStat) <- SummStatNames

  SummStat <- plyr::rbind.fill(SummStat, testStat)
  list("ArmData" = ArmData, "SummStat" = SummStat)
}

#------------------ -
# Compute the pair-wise test statistics for a given pair ID from summary Stat(df)
#------------------ -
pairTestStat <- function(pairID, df) {
  arms <- unlist(strsplit(pairID, split = "-"))
  SigmaEstimated <- T
  if (SigmaEstimated) {
    delta <- (df$Mean[df$Arm == arms[1]] - df$Mean[df$Arm == arms[2]])
    se.pair <- sqrt(df$SE[df$Arm == arms[1]]^2 / df$Completers[df$Arm == arms[1]] +
      df$SE[df$Arm == arms[2]]^2 / df$Completers[df$Arm == arms[2]])

    return(delta / se.pair)
  }
}

#----------------- -
# Generate Random Normal responses for the given (SimulationID,LookID,ArmID)
#----------------- -
normalResponse <- function(armID, n, mu, s, Arm.seed) {
  set.seed(Arm.seed)
  x <- rnorm(n, mu, s)
  x_bar <- mean(x)
  se.sample <- sd(x) * sqrt((n / (n - 1)))
  list("Arm" = armID, "Completers" = n, "Mean" = x_bar, "SE" = se.sample)
}

#----------------- -
# Summary Arm-Wise
#----------------- -
perArmData <- function(SimSeed, simID, lookID, mcpObj, gmcpSimObj) {
  if (lookID == 1) # for look 1 and FSD
    {
      SS.arm <- mcpObj$PlannedSS[1, ]
      mu.arm <- gmcpSimObj$Arms.Mean
      sd.arm <- gmcpSimObj$Arms.std.dev
      result_list <- lapply(1:length(SS.arm), function(x) {
        normalResponse(
          armID = names(SS.arm)[x], n = SS.arm[x],
          mu = mu.arm[x], s = sd.arm[x], Arm.seed = getRunSeed(SimSeed, simID, lookID, x)
        )
      })
      df <- data.frame(Arm = character(), Mean = numeric(), SE = numeric())
      df <- do.call(rbind, lapply(result_list, as.data.frame))
      rownames(df) <- NULL
      return(df)
    } else # interim and Final look
  {
    LookNrespose <- getInterimResposes(lookID, mcpObj, gmcpSimObj)
    mu.arm <- LookNrespose$mu.arm
    sd.arm <- LookNrespose$sd.arm
    SS.arm <- LookNrespose$SS.arm
    SS.arm <- SS.arm[names(mu.arm)]

    result_list <- lapply(1:length(SS.arm), function(x) {
      normalResponse(
        armID = names(SS.arm)[x], n = SS.arm[x],
        mu = mu.arm[x], s = sd.arm[x], Arm.seed = getRunSeed(SimSeed, simID, lookID, x)
      )
    })
    df <- data.frame(Arm = character(), Mean = numeric(), SE = numeric())
    df <- do.call(rbind, lapply(result_list, as.data.frame))
    rownames(df) <- NULL
    return(df)
  }
}

#---------------  -
# Response Generation for interim looks based on the available arms or re-allocated sample size(Implicit SSR)
#--------------- -
getInterimResposes <- function(lookID, mcpObj, gmcpSimObj) {
  planSS <- SS.arm <- mcpObj$PlannedSS[lookID, ] # If no Implicit SSR is needed SS.arm = Planned Sample Size

  if (gmcpSimObj$ImplicitSSR == "Selection") # Re-allocate samples from only de-selected arms
    {
      AditionalSS <- sum(planSS[which(mcpObj$ArmsRetained)])
      ss_frac <- gmcpSimObj$Arms.alloc.ratio[mcpObj$ArmsPresent] / sum(gmcpSimObj$Arms.alloc.ratio[mcpObj$ArmsPresent])
      SS.arm.add <- round(AditionalSS * ss_frac)
      SS.arm.add[1] <- AditionalSS - sum(SS.arm.add[-1])
      SS.arm[which(mcpObj$ArmsPresent)] <- SS.arm[which(mcpObj$ArmsPresent)] + SS.arm.add
      SS.arm[which(!mcpObj$ArmsPresent)] <- NA
    } else if (gmcpSimObj$ImplicitSSR == "All") # Allocate all the samples planed to the available arms
    {
      TotalAllocSS <- sum(planSS)
      ss_frac <- gmcpSimObj$Arms.alloc.ratio[mcpObj$ArmsPresent] / sum(gmcpSimObj$Arms.alloc.ratio[mcpObj$ArmsPresent])
      SS.arm[which(mcpObj$ArmsPresent)] <- round(TotalAllocSS * ss_frac)
      SS.arm[which(!mcpObj$ArmsPresent)] <- NA
    }

  mu.arm <- gmcpSimObj$Arms.Mean[mcpObj$ArmsPresent]
  sd.arm <- gmcpSimObj$Arms.std.dev[mcpObj$ArmsPresent]
  list("SS.arm" = SS.arm, "mu.arm" = mu.arm, "sd.arm" = sd.arm)
}

#-------------- -
# Generate seed for the given combination (Input Seed, SimulationID,LookID,ArmID)
getRunSeed <- function(SimSeed, simID, lookID, armIndex) {
  SimSeed + as.integer(paste(simID, lookID, armIndex, sep = ""))
}

#------------- -
# Compute Arm-Wise Look-Wise Panned Sample Size
#------------- -
getPlanSS <- function(gmcpSimObj) {
  ss_frac <- gmcpSimObj$Arms.alloc.ratio / sum(gmcpSimObj$Arms.alloc.ratio)
  ss_lk <- round(gmcpSimObj$Max_SS * gmcpSimObj$InfoFrac)
  ss_lk[gmcpSimObj$nLooks] <- gmcpSimObj$Max_SS
  ss_incr <- c(ss_lk[1], diff(ss_lk))
  ss_plan <- t(sapply(ss_incr, function(x) round(x * ss_frac)))

  if (gmcpSimObj$nLooks > 1) ss_plan[, 1] <- ss_incr - rowSums(ss_plan[, -1])

  ss_plan
}

#------------- -
# Identify retained arms
#------------- -
getDroppedArms <- function(ArmsPresent, rejflags, index_map) {
  SetR <- getArms(names(rejflags)[rejflags], index_map) # Arms falls under rejection group
  SetNR <- getArms(names(rejflags)[!rejflags], index_map) # Arms falls under non-rejection group
  SetD <- names(ArmsPresent)[!ArmsPresent] # Arms retained

  droppedArms <- union(setdiff(SetR, SetNR), SetD)
  ArmsPresent[droppedArms] <- F
  ArmsPresent
}

#------------- -
# Identify retained arms
#------------- -
getArmsPresent <- function(ArmsPresent, rejflags, HypoMap) {
  rejHypo <- HypoMap$Hypothesis[rejflags]
  nrejHypo <- HypoMap$Hypothesis[!rejflags]

  SetR <- getArms2(rejHypo, HypoMap) # Arms falls under rejection group
  SetNR <- getArms2(nrejHypo, HypoMap) # Arms falls under non-rejection group
  # SetD <- 1:length(ArmsPresent)[!ArmsPresent] #Arms retained

  # droppedArms <- union(setdiff(SetR,SetNR), SetD)
  ArmsPresent[setdiff(SetR, SetNR)] <- F
  ArmsPresent
}

#-------------- -
# Get Arms from Hypothesis index
#-------------- -
getArms2 <- function(SetH, HypoMap) {
  if (length(SetH) == 0) {
    return(NULL)
  } else {
    unique(unlist(HypoMap[
      HypoMap$Hypothesis %in% SetH,
      c("Control", "Treatment")
    ]))
  }
}

#-------------- -
# Get Arms from Hypothesis index
#-------------- -
getArms <- function(SetH, index_map) {
  if (length(SetH) == 0) {
    return(NULL)
  } else {
    return(unique(unlist(
      lapply(
        index_map$Names[index_map$GlobalIndexSet %in% SetH],
        function(x) {
          strsplit(x, split = "-")
        }
      )
    )))
  }
}

#-------------- -
# Inverse Normal Weights
#-------------- -
getInvNormWeights <- function(planSSIncr) {
  SS_looks_incr <- rowSums(planSSIncr)
  SS_looks_cumulative <- cumsum(SS_looks_incr)
  W_Norm <- matrix(NA, nrow = nrow(planSSIncr), ncol = nrow(planSSIncr))

  for (i in 1:nrow(W_Norm))
  {
    for (j in 1:i) W_Norm[i, j] <- sqrt(SS_looks_incr[j] / SS_looks_cumulative[i])
  }
  InvNormWeights <- W_Norm
  colnames(InvNormWeights) <- paste("W", 1:nrow(planSSIncr), sep = "")
  rownames(InvNormWeights) <- paste("Look", 1:nrow(planSSIncr), sep = "")
  W_Norm <- W_Norm[-1, ] # Removing the first row

  list(
    "InvNormWeightsTab" = knitr::kable(InvNormWeights, align = "c"),
    "W_Norm" = W_Norm
  )
}


#-------------- -
# boundaries for combining p-values method
#-------------- -
getPvalBdry <- function(alpha = 0.025,
                        nLooks = 3,
                        info_frac = c(1/3,2/3,1),
                        typeOfDesign = "asOF",
                        deltaWT = 0,
                        deltaPT1 = 0,
                        gammaA = 2,
                        userAlphaSpending = rpact::getDesignGroupSequential(
                          sided = 1, alpha = alpha,informationRates =info_frac,
                          typeOfDesign = "asOF")$alphaSpent)
  {
  UseExternal <- T
  if (UseExternal) # this part of the code can be replaced later with the internal R-codes
    {
    if(typeOfDesign == "WT"){
      des <- rpact::getDesignGroupSequential(
        kMax = nLooks, alpha = alpha,
        informationRates = info_frac,
        typeOfDesign = typeOfDesign,
        deltaWT = deltaWT
      )
    }else if(typeOfDesign == "PT"){
      des <- rpact::getDesignGroupSequential(
        kMax = nLooks, alpha = alpha,
        informationRates = info_frac,
        typeOfDesign = typeOfDesign,
        deltaPT1 = deltaPT1
      )
    }else if(typeOfDesign == "asHSD" || typeOfDesign == "asKD"){
      des <- rpact::getDesignGroupSequential(
        kMax = nLooks, alpha = alpha,
        informationRates = info_frac,
        typeOfDesign = typeOfDesign,
        gammaA = gammaA
      )
    }else if(typeOfDesign == "asUser"){
      des <- rpact::getDesignGroupSequential(
        kMax = nLooks, alpha = alpha,
        informationRates = info_frac,
        typeOfDesign = typeOfDesign,
        userAlphaSpending = userAlphaSpending
      )
    }else{
      des <- rpact::getDesignGroupSequential(
        kMax = nLooks, alpha = alpha,
        informationRates = info_frac,
        typeOfDesign = typeOfDesign
      )
    }

    Threshold <- des$stageLevels
    incr_alpha <- c(des$alphaSpent[1], diff(des$alphaSpent))

    # Boundary table
    bdryTab <- data.frame(
        "Look" = 1:nLooks,
        "Information_Fraction" = info_frac,
        "Incr_alpha_spent" = incr_alpha,
        "ZScale_Eff_Bbry" = des$criticalValues,
        "PValue_Eff_Bbry" = Threshold
      )

    colnames(bdryTab) <- c(
        "Looks", "InfoFrac", "Alpha(Incr.)",
        "Boundary(Z)", "Boundary(P-Value)"
      )
      bdryTab <- knitr::kable(bdryTab, align = "c")
    }
  list(
    "pValueBdryTab" = bdryTab,
    "Threshold" = Threshold
  )
}
