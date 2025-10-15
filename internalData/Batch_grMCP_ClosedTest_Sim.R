# File: Batch_grMCP_ClosedTest_Sim.R
# Code for calculating power of a MC design by performing closed testing using
# the raw p-values computed in each simulation

library(tidyverse)
library(graphicalMCP)
library(parallel)

# Function to run closed testing and build results dataframe
RunSimWithClosedTesting <- function(pvalues, graph, alpha, test_groups,
                                    n_hypotheses, test_type, n_eps, test_corr,
                                    bRunParallel = TRUE) {
  # Create column names for the results dataframe
  sAdjPColNames <- paste0("Adj_P_", 1:n_hypotheses)
  sRejColNames <- paste0("Rej_", 1:n_hypotheses)

  # Initialize the results dataframe
  dfResModel <- data.frame(matrix(NA, nrow = nrow(pvalues), ncol = 2*n_hypotheses + 2))

  # Set column names for the dataframe
  colnames(dfResModel) <- c(sAdjPColNames, sRejColNames, "RejAtLst1", "RejAll")

  # Helper function to process one simulation
  process_one_sim <- function(sim_index) {
    sim_results <- graphicalMCP::graph_test_closure(
      graph = graph,
      p = pvalues[sim_index, ],
      alpha = alpha,
      test_groups = test_groups,
      test_types = rep(test_type, n_eps),
      test_corr = rep(list(test_corr), n_eps)
    )

    # Adjust p-values for the hypotheses
    adj_p_vals <- sim_results$outputs$adjusted_p

    # Hypothesis rejection statuses
    rej_hypos <- sim_results$outputs$rejected * 1

    # Calculate RejAtLst1 and RejAll
    rej_at_lst1 <- ifelse(sum(rej_hypos) >= 1, 1, 0)
    rej_all <- prod(rej_hypos)

    # Return all results
    return(list(adj_p_vals = adj_p_vals,
                rej_hypos = rej_hypos,
                rej_at_lst1 = rej_at_lst1,
                rej_all = rej_all))
  }

  # If test_type is parametric, use parallel processing
  if (test_type == "parametric" && bRunParallel == TRUE) {
    cat("Using parallel processing for parametric test...\n")

    # Determine the number of cores to use (leave one for the OS)
    n_cores <- max(1, parallel::detectCores() - 1)

    # Create a cluster
    cl <- parallel::makeCluster(n_cores)

    # Export necessary objects to the cluster
    parallel::clusterExport(cl, varlist = c("graph", "alpha", "test_groups",
                                  "test_type", "n_eps", "test_corr", "pvalues"),
                  envir = environment())

    # Load required packages on each cluster node
    parallel::clusterEvalQ(cl, library(graphicalMCP))

    # Run in parallel
    all_results <- parallel::parLapply(cl, seq_len(nrow(pvalues)), process_one_sim)

    # Stop the cluster
    parallel::stopCluster(cl)

    # Fill the results dataframe
    for (nSim in seq_len(nrow(pvalues))) {
      result <- all_results[[nSim]]
      dfResModel[nSim, 1:n_hypotheses] <- result$adj_p_vals
      dfResModel[nSim, (n_hypotheses+1):(2*n_hypotheses)] <- result$rej_hypos
      dfResModel[nSim, "RejAtLst1"] <- result$rej_at_lst1
      dfResModel[nSim, "RejAll"] <- result$rej_all
    }
  } else {
    # Process sequentially for non-parametric tests
    for(nSim in seq_len(nrow(pvalues))) {
      # Use the same process_one_sim function for consistency
      result <- process_one_sim(nSim)

      # Store results in the dataframe
      dfResModel[nSim, 1:n_hypotheses] <- result$adj_p_vals
      dfResModel[nSim, (n_hypotheses+1):(2*n_hypotheses)] <- result$rej_hypos
      dfResModel[nSim, "RejAtLst1"] <- result$rej_at_lst1
      dfResModel[nSim, "RejAll"] <- result$rej_all
    }
  }

  return(dfResModel)
}

# Start timing the total execution
total_start_time <- Sys.time()
cat("Starting execution at:", format(total_start_time, "%Y-%m-%d %H:%M:%S"), "\n\n")

# Read raw p-values data
dfAllRawPValues <- read.csv("internalData/RawPValues_165_models.csv")

# Read graph input data
dfGraphInp <- read.csv("internalData/BatchInput_AGMCP_Sim_Bin_3.csv")

# We need to filter out models with Sidak test as graphicalMCP does not support it
dfSidak <- dfGraphInp %>% filter(test.type == "Sidak")

# Keeping only those raw p-value records that do not correspond to Sidak models
dfAllRawPValues <- dfAllRawPValues %>% filter(!(ModelID %in% dfSidak$ModelID))

# Similarly, keeping only those scenarios that do not correspond to Sidak models
dfGraphInp <- dfGraphInp %>% filter(test.type != "Sidak")

# Find unique ModelID values
uniqueModelIDs <- unique(dfAllRawPValues$ModelID)

# Split the dataframe by ModelID
dfRawPVals <- list()
for (id in uniqueModelIDs) {
  dfRawPVals[[as.character(id)]] <- dfAllRawPValues[dfAllRawPValues$ModelID == id, ]
}

# Filter rows corresponding to unique ModelIDs from dfAllRawPValues
dfGraphInp <- dfGraphInp[dfGraphInp$ModelID %in% uniqueModelIDs, ]

# Create arrays for WI and G for each ModelID
WI <- setNames(dfGraphInp$WI, dfGraphInp$ModelID)
G <- setNames(dfGraphInp$G, dfGraphInp$ModelID)

dAlphas <- dfGraphInp$alpha
nEPss <- dfGraphInp$nEps
nArmss <- dfGraphInp$nArms

# Flag for testing/debugging - set to TRUE for testing, FALSE for production
isTestMode <- FALSE  # Change this to FALSE for production use
nSims <- 0 # 100 # 10 # 1000
bParallel <- TRUE # Set to TRUE to run parametric test comps in parallel

# Flag for printing verbose output
bOutVerbose <- FALSE # Set to TRUE for verbose and FALSE for non-verbose output.

# Initialize list to store results
dfResults <- data.frame("ModelID" = uniqueModelIDs)
dfResults$PowAtLst1 <- rep(NA, length(uniqueModelIDs))
dfResults$RejectAll <- rep(NA, length(uniqueModelIDs))

nIters <- uniqueModelIDs # uniqueModelIDs[1]

# Loop through each unique ModelID
for (id in nIters) {
  # tStartTime <- Sys.time()
  print("-----------------------------------------------------------------")
  print("ModelID:")
  print(id)

  # Get the raw p-values for this ModelID
  fullPValues <- dfRawPVals[[as.character(id)]]

  # If in test mode, use only the first nSims rows, otherwise use all rows
  if (isTestMode) {
    pValues <- fullPValues[seq_len(min(nSims, nrow(fullPValues))), ]
    cat("TEST MODE: Using only the first", min(nSims, nrow(fullPValues)), "rows of p-values\n")
  } else {
    pValues <- fullPValues
  }

  pValues <- pValues %>% select(-c(ModelNum, ModelID))

  if(bOutVerbose) {
    print(paste("Num pvalues: ", nrow(pValues)))
    print(head(pValues))
  }

  # Number of endpoints and arms
  nEPs <- nEPss[dfGraphInp$ModelID==id]
  nArms <- nArmss[dfGraphInp$ModelID==id]
  nTrtms <- nArms - 1

  # Total and arm wise sample sizes
  nSampleSize <- dfGraphInp$SampleSize[dfGraphInp$ModelID==id]
  Arms.alloc.ratio <-
    eval(parse(text = dfGraphInp$Arms.alloc.ratio[dfGraphInp$ModelID==id]))
  dArmSS <- nSampleSize * Arms.alloc.ratio / sum(Arms.alloc.ratio)

  # Finding the distinct test groups
  lEpType <- eval(parse(text = dfGraphInp$lEpType[dfGraphInp$ModelID==id]))
  lTestGroups <- list()
  for (nEP in seq_along(lEpType)) {
    lTestGroups[[nEP]] <- (nArms-1) * (nEP-1) + (1:(nArms-1))
  }

  # alpha to be used for testing
  dAlpha <- dAlphas[dfGraphInp$ModelID==id]

  # Creating the graph object to be used for testing
  # Graph parameters
  weights <- eval(parse(text = WI[as.character(id)]))
  graphMatrix <- eval(parse(text = G[as.character(id)]))
  gGraph <- graph_create(hypotheses = weights, transitions = graphMatrix)

  if(bOutVerbose){
    print("weights:")
    print(weights)

    print("graphMatrix:")
    print(graphMatrix)

    print("graph:")
    print(gGraph)
  }

  sTest <- dfGraphInp$test.type[dfGraphInp$ModelID==id]

  # Mapping AdaptGMCP test type to graphicalMCP test type
  sAGMPTests <- c("Bonf", "Simes", "Dunnett")
  sgrMCPTests <- c("bonferroni", "simes", "parametric")

  sgrMCPTst <- sgrMCPTests[which(sAGMPTests == sTest)]

  if(bOutVerbose){
    print("test: ")
    print(sgrMCPTst)
  }

  # Correlation matrix for the test - corr.test
  mCorrTest <- if(sgrMCPTst == "parametric") {
    CalcCorrMatrixForTest(nEps = 1, nTrtms, dArmSS)
  } else {
    NA
  }

  if(bOutVerbose){
    print("corr.test: ")
    print(mCorrTest)
  }

  # Remove columns that contain all NAs
  na_columns <- sapply(pValues, function(x) all(is.na(x)))
  if(any(na_columns)) {
    if(bOutVerbose){
      cat("Removing", sum(na_columns), "columns with all NA values\n")
    }
    pValues <- pValues[, !na_columns, drop = FALSE]
  }

  pvalues <- as.matrix(pValues)

  # Use the function to compute the results dataframe
  dfResModel <- RunSimWithClosedTesting(pvalues = pvalues,
    graph = gGraph, alpha = dAlpha, test_groups = lTestGroups,
    n_hypotheses = nTrtms * nEPs, test_type = sgrMCPTst, n_eps = nEPs,
    test_corr = mCorrTest, bRunParallel = bParallel)

  # Print a summary of the results
  if(bOutVerbose){
    cat("Results summary:\n")
    cat("Number of simulations:", nrow(dfResModel), "\n")
    cat("Proportion with at least one rejection:", mean(dfResModel$RejAtLst1), "\n")
    cat("Proportion with all hypotheses rejected:", mean(dfResModel$RejAll), "\n")
  }

  # Store overall power results in dfResults
  dfResults[dfResults$ModelID==id, "PowAtLst1"] <- mean(dfResModel$RejAtLst1)
  dfResults[dfResults$ModelID==id, "RejectAll"] <- mean(dfResModel$RejAll)

  # # Save detailed results for this model ID
  # dfResModel$ModelID <- id
  #
  # # Create timestamp for output filename
  # timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  #
  # # Save results for this model to a CSV file
  # output_file <- paste0("internalData/DetailedResults_ModelID_", id, "_", timestamp, ".csv")
  # write.csv(dfResModel, output_file, row.names = FALSE)
  # cat("Detailed results saved to:", output_file, "\n")
  print("\n------------------------------------------------\n")
}

print(dfResults)
write_csv(dfResults, "internalData/grMCPPValSimOut.csv")

# Calculate and print the total execution time
total_end_time <- Sys.time()
total_execution_time <- difftime(total_end_time, total_start_time, units = "auto")

cat("\n--------------------------------------\n")
cat("Execution started at:", format(total_start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Execution ended at:  ", format(total_end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total execution time:", total_execution_time, attr(total_execution_time, "units"), "\n")
cat("--------------------------------------\n")

