# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

# The file contains supporting functions(Adaptations) for Adaptive GMCP Simulation Function AdaptGmcp_Simulation/adaptGMCP_SIM(.) ----
## Author: Ajoy.M

#------------ -
# Selection for a specific simulation
#------------ -
do_SelectionSim2 <- function(simID, mcpObj) {
  SelectFlag <- any(mcpObj$CurrentLook == mcpObj$SelectionLook)

  if (SelectFlag) {
    # Selection
    mcpObj <- getSelectedHypo2(simID=simID,
                               mcpObj = mcpObj)
  }
  mcpObj
}

#------------ -
# Identify Hypothesis based on selection rule
#------------- -
getSelectedHypo2 <- function(simID, mcpObj) {
  # The following script is applicable for Right Tail tests only#

  ############### Step-1: get the available hypothesis to select#################
  # Subset for the not rejected hypothesis
  # contHypo <- mcpObj$HypoMap[mcpObj$HypoPresent, ]
  # only exclude the treatments for which both endpoints are rejected
  hypothesis_not_rejected <- mcpObj$HypoMap[!mcpObj$rej_flag_Prev, ]
  contHypo <- mcpObj$HypoMap[mcpObj$HypoMap$Treatment %in% hypothesis_not_rejected$Treatment, ]
  # Selection based on end points
  if (mcpObj$SelectEndPoint != "overall") {
    contHypo <- contHypo[contHypo$Groups == mcpObj$SelectEndPoint, ]
  }
  #----------------------------------------------------------------------------

  ######## Step-2: get the scale variable on which selection will be performed##########
  if (mcpObj$SelectionScale == "pvalue") {
    ScaleVar <- mcpObj$SummStatDF[
      mcpObj$SummStatDF$LookID == mcpObj$CurrentLook,
      grep("RawPvalues", names(mcpObj$SummStatDF))
    ]
  } else if (mcpObj$SelectionScale == "delta") {
    ScaleVar <- mcpObj$SummStatDF[
      mcpObj$SummStatDF$LookID == mcpObj$CurrentLook,
      grep("Delta", names(mcpObj$SummStatDF))
    ]
  } else if (mcpObj$SelectionScale == "teststat") {
    ScaleVar <- mcpObj$SummStatDF[
      mcpObj$SummStatDF$LookID == mcpObj$CurrentLook,
      grep("TestStat", names(mcpObj$SummStatDF))
    ]
  } else if (mcpObj$SelectionScale == "stderror") {
    ScaleVar <- mcpObj$SummStatDF[
      mcpObj$SummStatDF$LookID == mcpObj$CurrentLook,
      grep("StdError", names(mcpObj$SummStatDF))
    ]
  }

  #-----------------------------------------------------------------------------
  ####################### Step-3: apply the selection rule ######################

  if(mcpObj$SelectionCriterion == "random"){
    ScaleVar <- NULL
  }else{
    ScaleVar <- as.numeric(ScaleVar[get_numeric_part(contHypo$Hypothesis)])

  }
  if (mcpObj$SelectionCriterion == "best") {
    # Best r
    if (mcpObj$SelectionScale == "stderror" || mcpObj$SelectionScale == "pvalue") {
      # for p-values and stderror the smaller the better
      ranks <- rank(ScaleVar, ties.method = "random")
    } else {
      # for test stat and delta the larger the better(right tail)
      ranks <- rank(-ScaleVar, ties.method = "random")
    }
    selectedH <- contHypo$Hypothesis[ranks <= mcpObj$SelectionParmeter]
  } else if (mcpObj$SelectionCriterion == "threshold") {
    # threshold based selection
    if (mcpObj$SelectionScale == "stderror" || mcpObj$SelectionScale == "pvalue") {
      # for p-values and stderror the smaller the better
      selectedH <- contHypo$Hypothesis[ScaleVar <= mcpObj$SelectionParmeter]
    } else {
      # for test stat and delta the larger the better(right tail)
      selectedH <- contHypo$Hypothesis[ScaleVar >= mcpObj$SelectionParmeter]
    }
  } else if (mcpObj$SelectionCriterion == "epsilon") {
    # epsilon neighbour of best
    if (mcpObj$SelectionScale == "stderror" || mcpObj$SelectionScale == "pvalue") {
      # for p-values and stderror the smaller the better
      ranks <- rank(ScaleVar, ties.method = "random")
    } else {
      # for test stat and delta the larger the better(right tail)
      ranks <- rank(-ScaleVar, ties.method = "random")
    }
    bestScale <- ScaleVar[which(ranks == 1)]
    selectedH <- contHypo$Hypothesis[ScaleVar <= bestScale + mcpObj$SelectionParmeter &
                                       ScaleVar >= bestScale - mcpObj$SelectionParmeter]

  } else if (mcpObj$SelectionCriterion == "random"){
    # Random Selection
    Seed <- getRunSeed(SimSeed = mcpObj$SimSeed,
                       simID = simID,
                       lookID = mcpObj$SelectionLook,
                       armIndex = 0)
    set.seed(seed = Seed)
    nSelecHyp <- sample(0:nrow(contHypo),size = 1)
    if(nSelecHyp == 0){
      selectedH <- c()
    }else{
      selectedH <- sample(contHypo$Hypothesis, nSelecHyp)
    }
  }

  #--------------------------------------------------------------------------------

  if (length(selectedH) != 0) # If the selection set is non empty
  {
    if (mcpObj$KeepAssosiatedEps) { # Keep all the associated hypothesis
      # notRejSet <- mcpObj$HypoMap[mcpObj$HypoPresent, ]
      notRejSet <- mcpObj$HypoMap[mcpObj$HypoMap$Treatment %in% hypothesis_not_rejected$Treatment, ]
      selectedArms <- getArms2(SetH = selectedH, HypoMap = notRejSet)
      AssociatedHypo <- notRejSet$Hypothesis[notRejSet$Treatment %in% selectedArms]

      selectedH <- AssociatedHypo
    }


    dropedArms <- setdiff(
      (1:length(mcpObj$ArmsPresent))[mcpObj$ArmsPresent],
      getArms2(SetH = selectedH, HypoMap = mcpObj$HypoMap)
    )

    mcpObj$ArmsRetained[dropedArms] <- T
    mcpObj$ArmsPresent[dropedArms] <- F
    mcpObj$SelectedIndex <- selectedH
    HypoPresentAfterSelct <- rep(F, length(mcpObj$HypoPresent))
    HypoPresentAfterSelct[get_numeric_part(selectedH)] <- T
    mcpObj$HypoPresent <- HypoPresentAfterSelct
    SelectedIDX <- rep(F, length(mcpObj$HypoPresent))
    SelectedIDX[get_numeric_part(mcpObj$SelectedIndex)] <- T
    mcpObj$DropedFlag <- (!SelectedIDX) & (!mcpObj$rej_flag_Curr) # Hypothesis not selected and not rejected are defined as droped

    mcpObj$IndexSet <- mcpObj$HypoMap$Hypothesis[mcpObj$HypoPresent]
    mcpObj$WH <- mcpObj$WH[which(apply(mcpObj$WH[mcpObj$IndexSet], 1, sum, na.rm = T) != 0), ]
  } else # If the selection set is empty
  {
    mcpObj$ContTrial <- F
  }
  mcpObj
}
