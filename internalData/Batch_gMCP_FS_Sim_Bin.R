# Batch_gMCP_FS_Sim_Bin.R
# File for executing a batch of gMCP simulation scenarios

library(tidyverse)
library(gMCP)
library(parallel)
source("internalData/gMCP_calcPower_wrapper.R", keep.source = TRUE)

# Start timing the total execution
total_start_time <- Sys.time()
cat("Starting execution at:", format(total_start_time, "%Y-%m-%d %H:%M:%S"), "\n\n")

###############################################################################
# CODE FOR CALCULATING POWER USING gMCP::graphTest() FUNCTION
# Read the entire data from the file "internalData/RawPValues_HPC.csv" into a dataframe dfAllRawPValues.
# Find unique values from the numeric column ModelID of this dataframe. Suppose there are "k" unique values.
# Split dfAllRawPValues by these unique values and store them in the array dfRawPVals of size k,
# where dfRawPVals[i] contains the data for ModelID = i.
# Now read the entire data from the file "internalData/BatchInput_FS_GMCP_Sim_Bin_2.csv" into a dataframe dfGraphInp.
# This dataframe also has a numeric ModelID column and string columns WI and G among other columns.
# Filter the rows from dfGraphInp corresponding to the unique ModelID values from dfAllRawPValues and create arrays of WI and G
# from these rows.
# Now for each unique ModelID = i, we have the dataframe dfRawPVals[i] and the scalar strings WI[i] and G[i].
# Create a loop going over the unique ModelID values.

# Create a graph object using WI and G cols for each unique ModelID
# For ModelID = i, Using dfRawPVals[i] and the ith graph, call gMCP::graphTest()
# Pass its output to gMCP::extractPower()

# Read raw p-values data
# dfAllRawPValues <- read.csv("internalData/RawPValues_HPC.csv")
dfAllRawPValues <- read.csv("internalData/RawPValues_165_models.csv")

# Find unique ModelID values
uniqueModelIDs <- unique(dfAllRawPValues$ModelID)

# Split the dataframe by ModelID
dfRawPVals <- list()
for (id in uniqueModelIDs) {
  dfRawPVals[[as.character(id)]] <- dfAllRawPValues[dfAllRawPValues$ModelID == id, ]
}

# Read graph input data
# dfGraphInp <- read.csv("internalData/BatchInput_AGMCP_Sim_Bin_2.csv")
dfGraphInp <- read.csv("internalData/BatchInput_AGMCP_Sim_Bin_3.csv")

# Filter rows corresponding to unique ModelIDs from dfAllRawPValues
dfGraphInp <- dfGraphInp[dfGraphInp$ModelID %in% uniqueModelIDs, ]

# Create arrays for WI and G for each ModelID
WI <- setNames(dfGraphInp$WI, dfGraphInp$ModelID)
G <- setNames(dfGraphInp$G, dfGraphInp$ModelID)

dAlphas <- dfGraphInp$alpha

# Flag for testing/debugging - set to TRUE for testing, FALSE for production
isTestMode <- FALSE # TRUE  # Change this to FALSE for production use
nSims <- 0 # 10 # 1000

# Flag for running parametric test computations in parallel
isParallelMode <- TRUE

# Initialize list to store results
dfResults <- data.frame("ModelID" = uniqueModelIDs)
dfResults$PowAtlst1 <- rep(NA, length(uniqueModelIDs))
dfResults$RejectAll <- rep(NA, length(uniqueModelIDs))

nIters <- uniqueModelIDs # uniqueModelIDs[1] #

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

  print(paste("Num pvalues: ", nrow(pValues)))
  print(head(pValues))

  # Evaluate strings to get actual values
  weights <- eval(parse(text = WI[as.character(id)]))
  graphMatrix <- eval(parse(text = G[as.character(id)]))
  dAlpha <- dAlphas[dfGraphInp$ModelID==id]

  # Create graph object
  graph <- new("graphMCP", m=graphMatrix, weights=weights)

  print("weights:")
  print(weights)

  print("graphMatrix:")
  print(graphMatrix)

  print("graph:")
  print(graph)

  sTest <- dfGraphInp[dfGraphInp$ModelID==id, ]$test.type

  print("test: ")
  print(sTest)

  # Correlation matrix for the test - corr.test
  mCorrTest <- if(sTest == "Dunnett" || sTest == "Partly-Parametric") {
    nEps <- dfGraphInp[dfGraphInp$ModelID==id, ]$nEps
    nTrtms <- dfGraphInp[dfGraphInp$ModelID==id, ]$nArms - 1
    nTotSS <- dfGraphInp[dfGraphInp$ModelID==id, ]$SampleSize
    nAllocRatio <- eval(parse(text =
                        dfGraphInp[dfGraphInp$ModelID==id, ]$Arms.alloc.ratio))
    nArmSS <- nTotSS * nAllocRatio / sum(nAllocRatio)

    CalcCorrMatrixForTest(nEps, nTrtms, nArmSS)
  } else {
    NULL
  }

  print("corr.test: ")
  print(mCorrTest)

  # Remove columns that contain all NAs
  na_columns <- sapply(pValues, function(x) all(is.na(x)))
  if(any(na_columns)) {
    cat("Removing", sum(na_columns), "columns with all NA values\n")
    pValues <- pValues[, !na_columns, drop = FALSE]
  }

  pvalues <- as.matrix(pValues)

  testResults <- if(isParallelMode) {
    if(sTest == "Dunnett" || sTest == "Partly-Parametric") {
      # Run graphTest in parallel and get combined results
      runParallelGraphTest(pvalues, graph, dAlpha, mCorrTest)
    } else {
      graphTest(pvalues, alpha = dAlpha, graph = graph, cr = mCorrTest)
    }
  } else {
    graphTest(pvalues, alpha = dAlpha, graph = graph,
                             cr = mCorrTest)
  }

  # Extract power from results
  out <- extractPower(testResults)  # Store results
  dfResults[dfResults$ModelID==id, ]$PowAtlst1 <- out$PowAtlst1
  dfResults[dfResults$ModelID==id, ]$RejectAll <- out$RejectAll

  print("Test output:")
  print(out)

#  tElapTime <- Sys.time() - tStartTime
#
#  cat("Elapsed time: ", tElapTime, "\n\n")
}

print(dfResults)

# Calculate and print the total execution time
total_end_time <- Sys.time()
total_execution_time <- difftime(total_end_time, total_start_time, units = "auto")

cat("\n--------------------------------------\n")
cat("Execution started at:", format(total_start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Execution ended at:  ", format(total_end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total execution time:", total_execution_time, attr(total_execution_time, "units"), "\n")
cat("--------------------------------------\n")

