# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

# The file contains supporting functions for Adaptive GMCP Simulation Function AdaptGmcp_Simulation/adaptGMCP_SIM(.) ----
## Author: Ajoy.M

#------------------ -
# Replace row-specific NA values from one data frame(df1) with the values from another(df2)
#------------------ -
fillNa <- function(rowID, df1, df2) {
  commonCol <- intersect(names(df1), names(df2))
  df1[rowID, commonCol] <- df2[rowID, commonCol]
  df1
}

#-------------- -
# Generate seed for the given combination (Input Seed, SimulationID,LookID,ArmID)
getRunSeed <- function(SimSeed, simID, lookID, armIndex, simID_Stage2 = 0) {
  # SimSeed + as.integer(paste(simID, simID_Stage2, lookID, armIndex, sep = ""))
  SimSeed + simID + simID_Stage2*10^5 + lookID *10 + armIndex
  # add an error handling condition for when this number becomes too high
  # and can't be used as an integer anymore
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
                          typeOfDesign = "asOF")$alphaSpent) {
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
