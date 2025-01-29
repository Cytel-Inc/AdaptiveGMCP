# The file contains supporting functions for Adaptive GMCP Simulation ----
## Author: Martin.Posh

#------------------ -
# Perform iterative simulation for Multi-Arm designs
#------------------ -

## Modification of AdaptGMCP to also work on MacOs/Linux

# Modify the function
modified_MAMSMEP_sim2 <- function (gmcpSimObj)
{
  preSimObjs <<- getPreSimObjs(gmcpSimObj = gmcpSimObj)
  starttime <- Sys.time()
  SummaryStatFile <- ArmWiseSummary <- PowerTab <- data.frame()
  powersName <- c("simID", "nG", "nC", "nD", "nF")
  PowerTab <- data.frame(matrix(nrow = 0, ncol = length(powersName)))
  EfficacyTable <- data.frame()
  names(PowerTab) <- powersName
  SelectionTab <- data.frame()
  if (gmcpSimObj$Method == "CER") {
    if (gmcpSimObj$Parallel) {
      cores <- parallel::detectCores()
      if(.Platform$OS.type == "windows")# check OS
        {
        cl <- parallel::makeCluster(cores[1] - 1, type = "PSOCK")
      parallel::clusterExport(cl, c("gmcpSimObj", "preSimObjs"))
      out <- parallel::parLapply(cl = cl, 1:gmcpSimObj$nSimulation, function(x){
                                  out_SingleSim <- SingleSimCER2(x, gmcpSimObj, preSimObjs)
                                  return(out_SingleSim)})
      parallel::stopCluster(cl)
      }
      else
        {
          ###--- code to run in parallel on macos
          out <- parallel::mclapply(1:gmcpSimObj$nSimulation, function(x) {
        out_SingleSim <- SingleSimCER2(x, gmcpSimObj, preSimObjs)
        return(out_SingleSim)
      },mc.cores = cores[1] - 1)
         ###---
    }
      }
    else {
      out <- lapply(1:gmcpSimObj$nSimulation, function(x) {
        SingleSimCER2(x, gmcpSimObj, preSimObjs)
      })
    }
  } else {
    if (gmcpSimObj$Parallel) {
      cores <- parallel::detectCores()
      if(.Platform$OS.type == "windows")#check OS
      {
        cl <- parallel::makeCluster(cores[1] - 1, type = "PSOCK")
        parallel::clusterExport(cl, c("gmcpSimObj", "preSimObjs"))
        out <- parallel::parLapply(cl = cl, 1:gmcpSimObj$nSimulation,
                                   function(x) {
                                     out_SingleSim <- SingleSimCombPValue2(x, gmcpSimObj,
                                                                           preSimObjs)
                                     return(out_SingleSim)
                                   })
      }
      else
      {
      ###--- code to run in parallel on macos
      out <- parallel::mclapply(1:gmcpSimObj$nSimulation, function(x) {
        out_SingleSim <- SingleSimCombPValue2(x, gmcpSimObj, preSimObjs)
        return(out_SingleSim)
      },mc.cores = cores[1] - 1)
      ###---
    }
      } else {
      out <- lapply(1:gmcpSimObj$nSimulation, function(x) {
        SingleSimCombPValue2(x, gmcpSimObj, preSimObjs)
      })
    }
  }

  # Aggregating results
  SuccessedSims <- 0
  for (i in 1:length(out)) {
    if (length(out[[i]]) == 1) {
      if (grepl(pattern = "Error", x = out[[i]])) {
        sprintf("Error Simulation %d ", i)
      }
    } else {
      SummaryStatFile <- rbind(SummaryStatFile, out[[i]]$SummStatDF)
      ArmWiseSummary <- rbind(ArmWiseSummary, out[[i]]$ArmWiseDF)
      PowerTab <- rbind(PowerTab, out[[i]]$powerCountDF)
      EfficacyTable <- plyr::rbind.fill(EfficacyTable, out[[i]]$EfficacyTable)
      SelectionTab <- rbind(SelectionTab, out[[i]]$SelectionDF)
      SuccessedSims <- SuccessedSims + 1
    }
  }

  # Processing simulation results
  Sim_power <- SimPowers(nSimulation = gmcpSimObj$nSimulation,
                         PowerTab = PowerTab)
  Sim_power_df <- Sim_power
  Sim_power <- knitr::kable(Sim_power,
                            align = "c",
                            col.names = c(' ','Overall Powers','Confidence Interval(95%)'))
  eff_count <- colSums(EfficacyTable[, -1])
  EffTab <- data.frame(Hypothesis = names(eff_count), Count = eff_count, Percentage = 100 * (eff_count / SuccessedSims), row.names = NULL)
  rownames(EffTab) <- NULL
  EffTab <- knitr::kable(EffTab, align = "c")

  if (gmcpSimObj$Selection) {
    SelcCount <- table(SelectionTab$SelectedHypothesis)
    SelcPerc <- 100 * (SelcCount / SuccessedSims)
    SelecTab <- data.frame(Hypothesis = names(SelcCount), Count = as.vector(SelcCount), Percentage = as.vector(SelcPerc), row.names = NULL)
    SelecTab <- knitr::kable(SelecTab, align = "c")
  } else {
    SelecTab <- NA
  }

  elapsedTime <- Sys.time() - starttime

  # Detailed output preparation
  detailOutput <- if (gmcpSimObj$Method == "CER") {
    ifElse(gmcpSimObj$SummaryStat,
           list(
             PlanSampleSizeCum = preSimObjs$planSS$CumulativeSamples,
             PlannedSigma = preSimObjs$Sigma,
             Boundary_Table = preSimObjs$plan_Bdry$PlanBdryTab,
             Overall_Powers = Sim_power,
             Overall_Powers_df = Sim_power_df,
             EfficacyTable = EffTab,
             SelectionTable = SelecTab,
             Summary_Stat = SummaryStatFile,
             ArmWiseSummary = ArmWiseSummary,
             Seed = preSimObjs$SimSeed,
             SuccessedSims = SuccessedSims,
             elapsedTime = elapsedTime
           ),
           list(
             PlanSampleSizeCum = preSimObjs$planSS$CumulativeSamples,
             PlannedSigma = preSimObjs$Sigma,
             Boundary_Table = preSimObjs$plan_Bdry$PlanBdryTab,
             Overall_Powers = Sim_power,
             Overall_Powers_df = Sim_power_df,
             EfficacyTable = EffTab,
             SelectionTable = SelecTab,
             Seed = preSimObjs$SimSeed,
             SuccessedSims = SuccessedSims,
             elapsedTime = elapsedTime
           )
    )
  } else {
    ifElse(gmcpSimObj$SummaryStat,
           list(
             PlanSampleSizeIncr = preSimObjs$planSS$IncrementalSamples,
             PlannedCorrelation = preSimObjs$PlanCorrelation,
             Boundary_Table = preSimObjs$pValBdry$pValueBdryTab,
             Inverse_Normal_Weights = preSimObjs$InvNormWeights$InvNormWeightsTab,
             Overall_Powers = Sim_power,
             Overall_Powers_df = Sim_power_df,
             EfficacyTable = EffTab,
             SelectionTable = SelecTab,
             Summary_Stat = SummaryStatFile,
             ArmWiseSummary = ArmWiseSummary,
             Seed = preSimObjs$SimSeed,
             SuccessedSims = SuccessedSims,
             elapsedTime = elapsedTime
           ),
           list(
             PlanSampleSize = preSimObjs$planSS$IncrementalSamples,
             PlannedCorrelation = preSimObjs$PlanCorrelation,
             Boundary_Table = preSimObjs$pValBdry$pValueBdryTab,
             Inverse_Normal_Weights = preSimObjs$InvNormWeights$InvNormWeightsTab,
             Overall_Powers = Sim_power,
             Overall_Powers_df = Sim_power_df,
             EfficacyTable = EffTab,
             SelectionTable = SelecTab,
             Seed = preSimObjs$SimSeed,
             SuccessedSims = SuccessedSims,
             elapsedTime = elapsedTime
           )
    )
  }

  # Return detailed output
  if (gmcpSimObj$plotGraphs) {
    list(DetailOutTabs = detailOutput, iniGraph = preSimObjs$iniGraph)
  } else {
    list(DetailOutTabs = detailOutput)
  }
}
