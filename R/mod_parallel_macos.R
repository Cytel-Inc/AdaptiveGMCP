# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

# The file contains supporting functions for Adaptive GMCP Simulation ----
## Author: Martin.Posh

#------------------ -
# Perform iterative simulation for Multi-Arm designs
#------------------ -

## Modification of AdaptGMCP to also work on MacOs/Linux

# Modify the function
#' @importFrom data.table data.table setDT := .N .SD setorder
#' @importFrom dplyr %>%
modified_MAMSMEP_sim2 <- function (gmcpSimObj)
{
  preSimObjs <<- getPreSimObjs(gmcpSimObj = gmcpSimObj)
  starttime <- Sys.time()
  SummaryStatFile <- ArmWiseSummary <- SelectionTab <- data.table()
  powersName <- c("simID", "nG", "nC", "nD", "nF")
  PowerTab <- data.table(matrix(nrow = 0, ncol = length(powersName)))
  EfficacyTable <- data.frame()
  data.table::setnames(PowerTab, powersName)
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
  startTime_postproc <- Sys.time()
  # Aggregating results
  SuccessedSims <- 0
  # Initialize stagewise rejection counters
  stage1Rejections <- 0
  stage2Rejections <- 0
  # rawpvalues <- list()
  # lRejStatus <- list()
  for (i in 1:length(out)) {
    if (length(out[[i]]) == 1) {
      if (grepl(pattern = "Error", x = out[[i]])) {
        sprintf("Error Simulation %d ", i)
      }
    } else {
      # rawpvalues[[i]] = out[[i]]$rawpvalues
      # lRejStatus[[i]] = out[[i]]$SummStatDF
      df <- out[[i]]$SummStatDF

      # Skip if dataframe is empty or NULL
      if (is.null(df) || nrow(df) == 0) next

      # Check if any RejStatus columns have TRUE in the first row
      rejCols <- grep("^RejStatus", names(df), value = TRUE)

      if (any(df[1, rejCols], na.rm = TRUE)) {
        stage1Rejections <- stage1Rejections + 1
      } else if (nrow(df) >= 2) {
        # If first row has no TRUE values and there's a second row,
        # check if any RejStatus columns have TRUE in the second row
        if (any(df[2, rejCols], na.rm = TRUE)) {
          stage2Rejections <- stage2Rejections + 1
        }
      }
      # SummaryStatFile <- data.table::rbindlist(list(SummaryStatFile, data.table(out[[i]]$SummStatDF)), use.names = TRUE, fill = TRUE)
      # ArmWiseSummary <- data.table::rbindlist(list(ArmWiseSummary, data.table(out[[i]]$ArmWiseDF)), use.names = TRUE, fill = TRUE)
      PowerTab <- data.table::rbindlist(list(PowerTab, data.table(out[[i]]$powerCountDF)), use.names = TRUE, fill = TRUE)
      # SelectionTab <- data.table::rbindlist(list(SelectionTab, data.table(out[[i]]$SelectionDF)), use.names = TRUE, fill = TRUE)
      SuccessedSims <- SuccessedSims + 1

    }
  }
  # save rawpvalues as rds with timestamp in the name
  # timestamp <- format(Sys.time(), "%y%m%d_%H%M%S")
  # saveRDS(rawpvalues, paste0("Debug_rawpvalues_", timestamp, ".rds"))
  # saveRDS(lRejStatus, paste0("Debug_rejstatus_", timestamp, ".rds"))
  # save(rawpvalues, lRejStatus, file = paste0("Debug_", timestamp, ".rda"))

  dfStagewiseRejections <- data.frame("Count" = c(stage1Rejections, stage2Rejections),
                                      "Percentage" = c(stage1Rejections, stage2Rejections)/gmcpSimObj$nSimulation)
  # EfficacyTable <- data.table::rbindlist(lapply(out, function(x) {
  #   if (is.list(x) && !is.null(x$EfficacyTable)) { # Check if EfficacyTable exists and is not NULL
  #     df <- as.data.frame(matrix(x$EfficacyTable,
  #                                nrow = nrow(x$EfficacyTable),
  #                                dimnames = dimnames(x$EfficacyTable)))
  #     return(df)
  #   }
  # })) # Fill missing values in the final combined table
  # cols <- grep("^RejStatus", names(EfficacyTable), value = TRUE)  # Identify relevant columns
  # EfficacyTable[
  #   , Hypothesis := {
  #     # Create a matrix of values for all rows and columns
  #     result_matrix <- do.call(cbind, lapply(seq_along(cols), function(i) {
  #       ifelse(.SD[[i]], sub("RejStatus", "H", cols[i]), "")
  #     }))
  #
  #     # Collapse each row into a single string
  #     apply(result_matrix, 1, function(row) {
  #       paste(Filter(nzchar, row), collapse = ", ")
  #     })
  #   },
  #   .SDcols = cols
  # ]
  # EfficacyTable <- EfficacyTable[
  #   Hypothesis != "",  # Exclude rows with blank Hypothesis
  #   .(Count = .N),     # Count rows per group
  #   by = Hypothesis    # Group by Hypothesis
  # ]
  #
  # EfficacyTable_totalRows <- sum(EfficacyTable$Count)
  # EfficacyTable[, Percentage := (Count / EfficacyTable_totalRows) * 100]
  # # Sort by Count in descending order
  # data.table::setorder(EfficacyTable, -Count)
  # Processing simulation results
  Sim_power <- SimPowers(nSimulation = gmcpSimObj$nSimulation,
                         nSimulation_Stage2 = gmcpSimObj$nSimulation_Stage2,
                         PowerTab = PowerTab)
  Sim_power_df <- Sim_power
  Sim_power <- knitr::kable(Sim_power,
                            align = "c",
                            col.names = c(' ','Overall Powers','Confidence Interval(95%)'))
  # eff_count <- colSums(EfficacyTable)
  # EffTab <- data.frame(Hypothesis = names(eff_count), Count = eff_count, Percentage = 100 * (eff_count / SuccessedSims), row.names = NULL)
  # rownames(EffTab) <- NULL
  # EffTab <- knitr::kable(EffTab, align = "c")
  # EffTab <- EfficacyTable

  # if (gmcpSimObj$Selection) {
  #   SelcCount <- table(SelectionTab$SelectedHypothesis)
  #   SelcPerc <- 100 * (SelcCount / SuccessedSims)
  #   SelecTab <- data.frame(Hypothesis = names(SelcCount), Count = as.vector(SelcCount), Percentage = as.vector(SelcPerc), row.names = NULL)
  #   SelecTab <- knitr::kable(SelecTab, align = "c")
  # } else {
  #   SelecTab <- NA
  # }

  elapsedTime <- Sys.time() - starttime
  elapsedTime_postProc <- Sys.time() - startTime_postproc
  gmcpSimObj$SummaryStat <- FALSE
  # Detailed output preparation
  detailOutput <- if (gmcpSimObj$Method == "CER") {
    ifElse(gmcpSimObj$SummaryStat,
           list(
             PlanSampleSizeCum = preSimObjs$planSS$CumulativeSamples,
             PlannedSigma = preSimObjs$Sigma,
             Boundary_Table = preSimObjs$plan_Bdry$PlanBdryTab,
             Overall_Powers = Sim_power,
             Overall_Powers_df = Sim_power_df,
             # EfficacyTable = EffTab,
             # SelectionTable = SelecTab,
             # Summary_Stat = SummaryStatFile,
             # ArmWiseSummary = ArmWiseSummary,
             Seed = preSimObjs$SimSeed,
             SuccessedSims = SuccessedSims,
             elapsedTime = elapsedTime,
             elapsedTime_postProc = elapsedTime_postProc,
             power_raw = PowerTab,
             stagewiseRejections = dfStagewiseRejections
           ),
           list(
             PlanSampleSizeCum = preSimObjs$planSS$CumulativeSamples,
             PlannedSigma = preSimObjs$Sigma,
             Boundary_Table = preSimObjs$plan_Bdry$PlanBdryTab,
             Overall_Powers = Sim_power,
             Overall_Powers_df = Sim_power_df,
             # EfficacyTable = EffTab,
             # SelectionTable = SelecTab,
             Seed = preSimObjs$SimSeed,
             SuccessedSims = SuccessedSims,
             elapsedTime = elapsedTime,
             elapsedTime_postProc = elapsedTime_postProc,
             power_raw = PowerTab,
             stagewiseRejections = dfStagewiseRejections
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
             # EfficacyTable = EffTab,
             # SelectionTable = SelecTab,
             # Summary_Stat = SummaryStatFile,
             # ArmWiseSummary = ArmWiseSummary,
             Seed = preSimObjs$SimSeed,
             SuccessedSims = SuccessedSims,
             elapsedTime = elapsedTime,
             power_raw = PowerTab,
             stagewiseRejections = dfStagewiseRejections
           ),
           list(
             PlanSampleSize = preSimObjs$planSS$IncrementalSamples,
             PlannedCorrelation = preSimObjs$PlanCorrelation,
             Boundary_Table = preSimObjs$pValBdry$pValueBdryTab,
             Inverse_Normal_Weights = preSimObjs$InvNormWeights$InvNormWeightsTab,
             Overall_Powers = Sim_power,
             Overall_Powers_df = Sim_power_df,
             # EfficacyTable = EffTab,
             # SelectionTable = SelecTab,
             Seed = preSimObjs$SimSeed,
             SuccessedSims = SuccessedSims,
             elapsedTime = elapsedTime,
             power_raw = PowerTab,
             stagewiseRejections = dfStagewiseRejections
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
