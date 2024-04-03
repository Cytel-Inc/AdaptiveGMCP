# The file contains supporting functions(Adaptations) for Adaptive GMCP Simulation Function AdaptGmcp_Simulation/adaptGMCP_SIM(.) ----
## Author: Ajoy.M

#------------ -
# Selection for a specific simulation
#------------ -
do_SelectionSim <- function(mcpObj, gmcpSimObj, PreSimObj) {
  SelectFlag <- any(mcpObj$CurrentLook == gmcpSimObj$SelectionLook)

  if (SelectFlag) {
    mcpObj$SelectionLook <- c(mcpObj$SelectionLook, (mcpObj$CurrentLook + 1))

    # Selection
    mcpObj <- getSelectedHypo(mcpObj = mcpObj, gmcpSimObj = gmcpSimObj, PreSimObj = PreSimObj)

    for (i in 1:length(mcpObj$DropedFlag))
    {
      if (!mcpObj$DropedFlag[i] & !mcpObj$rej_flag_Prev[i]) {
        mcpObj$DropedFlag[i] <- !any(mcpObj$SelectedIndex == names(mcpObj$DropedFlag[i]))
      }
    }

    mcpObj$IndexSet <- intersect(mcpObj$IndexSet, mcpObj$SelectedIndex)
    mcpObj$WH <- mcpObj$WH[which(apply(mcpObj$WH[mcpObj$IndexSet], 1, sum, na.rm = T) != 0), ]
  }
  mcpObj
}

#------------ -
# Selection for a specific simulation
#------------ -
do_SelectionSim2 <- function(mcpObj) {
  SelectFlag <- any(mcpObj$CurrentLook == mcpObj$SelectionLook)

  if (SelectFlag) {
    # Selection
    mcpObj <- getSelectedHypo2(mcpObj = mcpObj)
  }
  mcpObj
}

#------------ -
# Identify Hypothesis based on selection rule
#------------- -
getSelectedHypo2 <- function(mcpObj) {
  # The following script is applicable for Right Tail tests only#

  ############### Step-1: get the available hypothesis to select#################
  # Subset for the not rejected hypothesis
  contHypo <- mcpObj$HypoMap[mcpObj$HypoPresent, ]

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
  ScaleVar <- as.numeric(ScaleVar[get_numeric_part(contHypo$Hypothesis)])

  #-----------------------------------------------------------------------------
  ####################### Step-3: apply the selection rule ######################
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
  }

  #--------------------------------------------------------------------------------

  if (length(selectedH) != 0) # If the selection set is non empty
    {
      if (mcpObj$KeepAssosiatedEps) { # Keep all the associated hypothesis
        notRejSet <- mcpObj$HypoMap[mcpObj$HypoPresent, ]
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


#------------ -
# Identify Hypothesis based on selection rule
#------------- -
getSelectedHypo <- function(mcpObj, gmcpSimObj, PreSimObj) {
  df1 <- PreSimObj$index_map[PreSimObj$index_map$GlobalIndexSet %in% mcpObj$IndexSet, ] # for the present hypothesis
  if (
    gmcpSimObj$SelectionScale == "pvalue" ## Selection based on hypothesis
  ) {
    # Remaining Hypothesis after rejection
    ScaleVar <- mcpObj$SummStat[, paste("RawPvalues", get_numeric_part(mcpObj$IndexSet), sep = "")]

    if (length(mcpObj$IndexSet) == 1) names(ScaleVar) <- paste("RawPvalues", get_numeric_part(mcpObj$IndexSet), sep = "")

    ScaleVal <- as.numeric(ScaleVar)
    names(ScaleVal) <- names(ScaleVar)

    if (gmcpSimObj$SelectionCriterion == "best") {
      # the least 1:gmcpSimObj$SelectionParmeter hypothesis is selected based on p-values
      SelecIndex <- get_numeric_part(names(sort(abs(ScaleVal), decreasing = F)[1:gmcpSimObj$SelectionParmeter]))
      if (length(SelecIndex) != 0) {
        selectedH <- paste("H", SelecIndex, sep = "")
      } else {
        selectedH <- c()
      }
    } else if (gmcpSimObj$SelectionCriterion == "threshold") {
      # the p-values less than the gmcpSimObj$SelectionParmeter threshold will be selected

      SelecIndex <- get_numeric_part(names(ScaleVal <= gmcpSimObj$SelectionParmeter))
      if (length(SelecIndex) != 0) {
        selectedH <- paste("H", SelecIndex, sep = "")
      } else {
        selectedH <- c()
      }
    }
    if (length(selectedH) != 0) # If the selection set is non empty
      {
        dropedArms <- setdiff(names(mcpObj$ArmsPresent[mcpObj$ArmsPresent]), getArms(selectedH, PreSimObj$index_map))
        mcpObj$ArmsRetained[dropedArms] <- T
        mcpObj$ArmsPresent[dropedArms] <- F
        mcpObj$SelectedIndex <- selectedH
      } else # If the selection set is empty
    {
      mcpObj$ContTrial <- F
    }
  } else if ( # Selection based on Arms
    gmcpSimObj$SelectionScale == "ArmID"
  ) {
    if (is.character(gmcpSimObj$SelectionParmeter)) {
      selectedArms <- gmcpSimObj$SelectionParmeter
    } else {
      selectedArms <- names(gmcpSimObj$ArmID[gmcpSimObj$ArmID %in% gmcpSimObj$SelectionParmeter])
    }

    dropedArms <- setdiff(names(mcpObj$ArmsPresent[mcpObj$ArmsPresent]), selectedArms)

    selectedH <- df1$GlobalIndexSet[unlist(lapply(df1$Names, function(x) {
      !any(unlist(strsplit(x, split = "-")) %in% dropedArms)
    }))]

    if (length(selectedH) != 0) # If the selection set is non empty
      {
        mcpObj$ArmsRetained[dropedArms] <- T
        mcpObj$ArmsPresent[dropedArms] <- F
        mcpObj$SelectedIndex <- selectedH
      } else # If the selection set is empty
    {
      mcpObj$ContTrial <- F
    }
  }
  mcpObj
}
