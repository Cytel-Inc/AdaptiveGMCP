# File: graphicalMCP_pow_calc_wrapper.R
# File containing a wrapper function for the function graph_calculate_power()
# from the graphicalMCP package and related utility functions

library(graphicalMCP)
library(purrr)
library(tidyverse)

# Function for simulating a batch of test cases
graphicalMCP_Wrapper <- function(InputDF) {
  # Update the dataframe column names in the following mapping in case
  # the names in the input csv/excel changes
  lOut <- list()

  for (nModelNum in 1:nrow(InputDF)) {
    # Start timer for this iteration
    start_time <- Sys.time()
    out <- tryCatch(
      {
        grphMCP_run1TestCase(InputDF = InputDF[nModelNum, ])
      },
      error = function(err) {
        paste0("Model ", nModelNum, " execution failed.")
      }
    )

    # End timer and calculate time taken in seconds
    end_time <- Sys.time()
    time_taken <- as.numeric(difftime(end_time, start_time, units = "hours"))

    OutTab <- data.frame("ModelID" = InputDF[nModelNum, "ModelID"])

    OutTab$ExpRejections <- out$power$rejection_expected
    OutTab$PowAtlst1 <- out$power$power_at_least_1
    OutTab$PowAll <- out$power$power_all

    # Add TimeTaken column
    OutTab$HoursTaken <- time_taken

    # add each iterations power table to a list
    lOut[[nModelNum]] <- OutTab

    passedTxt <- paste0("Model ", nModelNum, " execution completed successfully.")
    cat("\n", passedTxt, "\n")
    print(paste0("Power table for ", nModelNum, ":"))
    print(OutTab)
  }

  # rbind power tables for each iteration to produce a single table
  dfOut <- do.call(rbind, lOut)
  dfOut <- dplyr::left_join(dfOut, InputDF, by = "ModelID")

  return(dfOut)
}

# Function for simulating one test case
grphMCP_run1TestCase <- function(InputDF) {
  # mapping to link simMAMSMEP function arguments with csv columns
#  Method <- InputDF$Method
  nSampleSize <- InputDF$SampleSize
  dAlpha <- InputDF$alpha
  # TestStatCon <- InputDF$TestStatCon
  # TestStatBin <- InputDF$TestStatBin
  # FWERControl <- InputDF$FWERControl
  nArms <- InputDF$nArms
  nTrtms <- nArms-1
  nEps <- InputDF$nEps
  lEpType <- eval(parse(text =InputDF$lEpType))
  Arms.Mean <- eval(parse(text = InputDF$Arms.Mean))
  Arms.Prop <- eval(parse(text = InputDF$Arms.Prop))
  Arms.std.dev <- eval(parse(text = InputDF$Arms.std.dev))
  Arms.alloc.ratio <- eval(parse(text = InputDF$Arms.alloc.ratio))
  EP.Corr <- eval(parse(text = InputDF$EP.Corr))
  WI <- eval(parse(text = InputDF$WI))
  G <- eval(parse(text = InputDF$G))
  test.type <- InputDF$test.type
  # info_frac <- eval(parse(text = InputDF$info_frac))
  # typeOfDesign <- ifelse(is.na(InputDF$typeOfDesign), "asOF", InputDF$typeOfDesign)
  # MultipleWinners <- InputDF$MultipleWinners
  # MultipleWinners <- ifelse(is.na(MultipleWinners),FALSE,MultipleWinners)
  # Selection <- InputDF$Selection
  # Selection <- ifelse(is.na(Selection),FALSE, Selection)
  # CommonStdDev <- ifelse(is.na(InputDF$CommonStdDev), F, InputDF$CommonStdDev)
  # SelectionLook <- InputDF$SelectionLook
  # SelectEndPoint <- InputDF$SelectEndPoint
  # SelectionScale <- InputDF$SelectionScale
  # SelectionCriterion <- InputDF$SelectionCriterion
  # SelectionParmeter <- InputDF$SelectionParmeter
  # KeepAssosiatedEps <- InputDF$KeepAssosiatedEps
  # ImplicitSSR <- InputDF$ImplicitSSR
  # ImplicitSSR <- ifelse(is.na(ImplicitSSR),FALSE,ImplicitSSR)
  nSimulation <- InputDF$nSimulation
  Seed <- InputDF$Seed
  # SummaryStat <- InputDF$SummaryStat
  # plotGraphs <- InputDF$plotGraphs
  # Parallel <- InputDF$Parallel
  # nSimulation_Stage2 <- InputDF$nSimulation_Stage2
  # put the following code in try catch so the loop continues even if one iteration fails
  Seed <- if (!is.na(suppressWarnings(as.numeric(Seed)))) as.numeric(Seed) else Seed

  # Arm wise sample sizes
  dArmSS <- nSampleSize * Arms.alloc.ratio / sum(Arms.alloc.ratio)

  lMargPows <- list()
  lTestGroups <- list()

  for (nEP in 1:length(lEpType)) {
    sEP <- lEpType[[nEP]]

    lArmParams <- if(sEP == "Binary") {
      list(Props = Arms.Prop[[nEP]])
    } else if(sEP == "Continuous") {
      list(Means = Arms.Mean[[nEP]], StdDevs = Arms.std.dev[[nEP]])
    } else {
      stop("Invalid endpoint type encountered!")
    }

    # Calculating marginal powers for the hypo corresponding to this endpoint
    lMargPows[[nEP]] <- CalcMarginalPows(sEP, nArms, lArmParams, dArmSS, dAlpha)

    nHypos <- (nArms-1) * (nEP-1) + (1:(nArms-1))
    names(lMargPows[[nEP]]) <- paste0("H", nHypos)

    lTestGroups[[nEP]] <- nHypos
  }

  # Creating the graph object to be used for testing
  gGraph <- graph_create(hypotheses = WI, transitions = G)

  # Calculating the corr.sim matrix
  mCorrSim <- CalcCorrMatrixForSim(nEps, nTrtms, dArmSS, EP.Corr)

  # Correlation matrix for the test - corr.test
  mCorrTest <- if(test.type == "Dunnett" || test.type == "Partly-Parametric") {
  #  CalcCorrMatrixForTest(nEps, nTrtms, dArmSS)
    CalcCorrMatrixForTest(nEps = 1, nTrtms, dArmSS)
  } else {
    NA
  }

  # Flattening the marginal powers list as needed by graphicalMCP
  dMargPowVec <- unlist(lMargPows)

  # Mapping AdaptGMCP test type to graphicalMCP test type
  sAGMPTests <- c("Bonf", "Simes", "Dunnett", "Sidak")
  sgrMCPTests <- c("bonferroni", "simes", "parametric", "sidak")

  sgrMCPTst <- sgrMCPTests[which(sAGMPTests == test.type)]

  lPow <- graphicalMCP::graph_calculate_power(
    graph = gGraph, alpha = dAlpha, power_marginal = dMargPowVec,
#    test_groups = list(seq_along(nEps*nTrtms)), test_types = sgrMCPTst,
    test_groups = lTestGroups, test_types = rep(sgrMCPTst, nEps),
    test_corr = rep(list(mCorrTest), nEps), sim_n = nSimulation,
    sim_corr = mCorrSim)

  return(lPow)
}

# 1. Function for calculating marginal powers
# 2. Function for calculating non-centrality parameters for binomial distribution
# 3. Function for calculating non-centrality parameters for normal distribution
# 4. Use CalcCorrMatrixForSim() from gMCP_calcPower_wrapper.R for calculating the
#    correlation matrix for simulating test stats.
# 5. Use graphicalMCP::graph_create() for creating a graph using the node weights
#    and transition matrix.

# Function for calculating marginal powers
# sEPType: endpoint type "Binary" or "Continuous"
# nArms: number of arms including control
# lArmParams: list specifying arm wise response parameters
#             In case of continuous endpoints, must contain a vector of means
#             and a vector of std dev (both arm wise).
#             In case of binary endpoints, must contain a vector of arm wise
#             proportions.
# dArmSS: vector of arm wise sample sizes
# dAlpha: alpha to be used for two-sample one-sided test (must be same as the
#         alpha used for the MCP)
CalcMarginalPows <- function(sEPType, nArms, lArmParams, dArmSS, dAlpha) {
  # Calculating non-centrality parameters based on endpoint type
  dNCParams <- if(sEPType == "Continuous") {
    CalcNonCentrParamForNorm(nArms, lArmParams$Means, lArmParams$StdDevs,
                             dArmSS, dAlpha)
  } else if (sEPType == "Binary") {
    CalcNonCentrParamForBin(nArms, lArmParams$Props, dArmSS, dAlpha)
  } else {
    stop(sEPType, ": Invalid endpoint type encountered!")
  }

  # Calculating the marginal powers now
  dMargPows <- pnorm(qnorm(dAlpha, lower.tail = F),
                     mean = dNCParams, sd = 1, lower.tail = F)

  return(dMargPows)
}

# Function for calculating non-centrality parameters for binomial distribution
# nArms: number of arms including control
# dArmProps: vector of arm wise proportions
# dArmSS: vector of arm wise sample sizes
# dAlpha: alpha to be used for two-sample one-sided test (must be same as the
#         alpha used for the MCP)
CalcNonCentrParamForBin <- function(nArms, dArmProps, dArmSS, dAlpha) {
  dUnpooledVar <- dArmProps[-1] * (1 - dArmProps[-1]) / dArmSS[-1] +
    dArmProps[1] * (1 - dArmProps[1]) / dArmSS[1]

  dNonCentPar <- (dArmProps[-1] - dArmProps[1]) / sqrt(dUnpooledVar)

  nSign <- map_int(dArmProps[-1],
                   function(p, p1) { if(p >= p1) {1} else {-1} },
                   dArmProps[1])

  return(nSign * dNonCentPar)
}

# Function for calculating non-centrality parameters for normal distribution
# nArms: number of arms including control
# dArmMeans: vector of arm wise means
# dArmSDs: vector of arm wise std dev
# dArmSS: vector of arm wise sample sizes
# dAlpha: alpha to be used for two-sample one-sided test (must be same as the
#         alpha used for the MCP)
CalcNonCentrParamForNorm <- function(nArms, dArmMeans, dArmSDs,
                                     dArmSS, dAlpha) {
  dVariance <- dArmSDs[-1]^2 / dArmSS[-1] + dArmSDs[1]^2 / dArmSS[1]
  dNonCentPar <- (dArmMeans[-1] - dArmMeans[1]) / sqrt(dVariance)

  nSign <- map_int(dArmMeans[-1],
                   function(m, m1) { if(m >= m1) {1} else {-1} },
                   dArmMeans[1])

  return(nSign * dNonCentPar)
}
