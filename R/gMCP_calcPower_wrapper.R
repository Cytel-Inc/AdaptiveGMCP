# gMCP_calcPower_wrapper.R
# File containing a wrapper function for the gMCP::calcPower() function so that
# it can be called in a batch

library(gMCP)
library(dplyr)
library(tidyr)
library(purrr)

#' Wrapper function that takes all inputs through a single R dataframe and
#' calls gMCP::calcPower() on each row of the dataframe (one test case)
#' @param InputDF R Dataframe: This is the csv/excel input data in the R dataframe format
gMCP_CalcPowerWrapper <- function(InputDF) {
  # Update the dataframe column names in the following mapping in case
  # the names in the input csv/excel changes
  lOut <- list()
  for (nModelNum in 1:nrow(InputDF)) {
    # Start timer for this iteration
    start_time <- Sys.time()
    out1 <- tryCatch(
      {
        gMCP_Run1TestCase(InputDF = InputDF[nModelNum, ])
      },
      error = function(err) {
        paste0("Model ", nModelNum, " execution failed.")
      }
    )

    out <- out1$gMCPOut
    inp <- out1$gMCPInp

    dfGMCPInp <- FlattenGMCPInp(inp)

    # End timer and calculate time taken in seconds
    end_time <- Sys.time()
    time_taken <- as.numeric(difftime(end_time, start_time, units = "hours"))

    errorLog <- c(!grepl("Invalid", out[[1]]), !is.character(out), !is.null(out))

    if (!any(errorLog == F)) {
      OutTab <- data.frame("ModelID" = InputDF[nModelNum, "ModelID"])

      OutTab$ExpRejections <- out$ExpRejections
      OutTab$PowAtlst1 <- out$PowAtlst1

      # locPow <- t(out$LocalPower)
      # colnames(locPow) <- paste0("LocalPower_", colnames(locPow))
      #
      # OutTab <- cbind(OutTab, locPow)

      OutTab$HoursTaken <- time_taken

      OutTab <- cbind(OutTab, dfGMCPInp)

      # add each iterations power table to a list
      lOut[[nModelNum]] <- OutTab

      passedTxt <- paste0("Model ", nModelNum, " execution completed successfully.")
      cat("\n", passedTxt, "\n")
      print(paste0("Power table for ", nModelNum, ":"))
      print(OutTab)
    } else if (grepl("Invalid", out[[1]])) {
      failTxt <- paste0("Model ", nModelNum, " execution failed.")
      cat("\n", failTxt, "\n")
      cat("\n Details \n")
      print(unlist(out))
    } else {
      print(out)
    }
  }

  # rbind power tables for each iteration to produce a single table
  dfOut <- do.call(rbind, lOut)
  dfOut <- dplyr::left_join(dfOut, InputDF, by = "ModelID")

  # dfOut <- cbind(dfOut, dfGMCPInp)

  # lOut <- list(gMCPOut = dfOut, gMCPInp = inp)
  # return(lOut)

  return(dfOut)
}


gMCP_Run1TestCase <- function(InputDF) {
  # Extracting input
  mGraphMat <- eval(parse(text = InputDF$G))
  dWeights <- eval(parse(text = InputDF$WI))
  dAlpha <- InputDF$alpha

#  dArmProps <- unlist(eval(parse(text = InputDF$Arms.Prop)))
  lArmProps <- eval(parse(text = InputDF$Arms.Prop))

  SampleSize <- InputDF$SampleSize
  Arms.alloc.ratio <- eval(parse(text = InputDF$Arms.alloc.ratio))
  nSim <- InputDF$nSimulation
  nArms <- InputDF$nArms
  nEps <- InputDF$nEps
  EP.Corr <- eval(parse(text = InputDF$EP.Corr))
  test.type <- InputDF$test.type

  # Deriving other input parameters from the specified input
  nTrtms <- nArms-1

  # MCP graph
  gGraph <- new("graphMCP", m=mGraphMat, weights=dWeights)

  # Arm wise sample sizes
  nArmSS <- SampleSize * Arms.alloc.ratio / sum(Arms.alloc.ratio)

  # Z score values for the binomial MCP problem
  # dZMeans <- (dArmProps[-1] - dArmProps[1]) /
  #   sqrt(dArmProps[-1]*(1-dArmProps[-1])/nArmSS[-1] +
  #          dArmProps[1]*(1-dArmProps[1])/nArmSS[1])

  lZMeans <- map(lArmProps, function(x, nArmSS) {
    return(
      (x[-1] - x[1]) / sqrt(x[-1]*(1-x[-1])/nArmSS[-1] + x[1]*(1-x[1])/nArmSS[1])
    )
  }, nArmSS)

  dZMeans <- unlist(lZMeans)

  # Correlation matrix for the test - corr.test
  mCorrTest <- if(test.type == "Dunnett" || test.type == "Partly-Parametric") {
    CalcCorrMatrixForTest(nEps, nTrtms, nArmSS)
  } else {
    NULL
  }

  # Calculating the corr.sim matrix
  mCorrSim <- CalcCorrMatrixForSim(nEps, nTrtms, nArmSS, EP.Corr)

  # Finally calling calcPower() - COMMENTED OUT
  # out <- calcPower(graph=gGraph, alpha=dAlpha, mean=dZMeans, corr.sim = mCorrSim,
  #                  corr.test = mCorrTest, n.sim = nSim)

  # Using our dummy function instead
  out <- DummyCalcPower(graph=gGraph, alpha=dAlpha, mean=dZMeans, corr.sim = mCorrSim,
                      corr.test = mCorrTest, n.sim = nSim)

  lInp <- list(graph=gGraph, alpha=dAlpha, mean=dZMeans, corr.sim = mCorrSim,
               corr.test = mCorrTest, n.sim = nSim)
  lOut <- list(gMCPOut = out, gMCPInp = lInp)

  return(lOut)
}

# Function for calculating the correlation matrix used for testing
CalcCorrMatrixForTest <- function(nEps, nTrtms, nArmSS){
  mCorrTest <- matrix(nrow = nEps*nTrtms, ncol = nEps*nTrtms)
  diag(mCorrTest) <- 1

  # Control arm sample size
  nCtrlSS <- nArmSS[1]

  for (m in 1:nEps) {
    for (i in 1:(nTrtms-1)) {
      for (j in (i+1):nTrtms){
        SS_i <- nArmSS[i+1]
        SS_j <- nArmSS[j+1]

        mCorrTest[(m-1)*nTrtms+i, (m-1)*nTrtms+j] <-
          mCorrTest[(m-1)*nTrtms+j, (m-1)*nTrtms+i] <-
          sqrt((SS_i/(SS_i+nCtrlSS))*(SS_j/(SS_j+nCtrlSS)))
      }
    }
  }

  return(mCorrTest)
}

# Function for calculating the correlation matrix used for simulating the Z
# statistics in gMCP
CalcCorrMatrixForSim <- function(nEps, nTrtms, nArmSS, EP.Corr) {
  mCorrSim <- matrix(nrow = nEps*nTrtms, ncol = nEps*nTrtms)
  diag(mCorrSim) <- 1

  # Control arm sample size
  nCtrlSS <- nArmSS[1]

  for (p in 1:(nEps*nTrtms-1)) {
    for (q in (p+1):(nEps*nTrtms)) {
      nRowEP <- ceiling(p/nTrtms)
      nRowTrt <- p - nTrtms*(nRowEP-1)

      nColEP <- ceiling(q/nTrtms)
      nColTrt <- q - nTrtms*(nColEP-1)

      nRowSS <- nArmSS[nRowTrt+1]
      nColSS <- nArmSS[nColTrt+1]

#      dEPCorr <- ifelse(is.na(EP.Corr), 0, EP.Corr[nRowEP, nColEP])

      dEPCorr <- if(is.matrix(EP.Corr)) {
        EP.Corr[nRowEP, nColEP]
      } else {
        0
      }

      if(nRowEP == nColEP){ # Same endpoint
        # Note: Since we are looking at strictly off-diagonal entries,
        # nRowTrt and nColTrt have to be unequal.
        mCorrSim[p, q] <- mCorrSim[q, p] <-
          sqrt((nRowSS/(nRowSS+nCtrlSS))*(nColSS/(nColSS+nCtrlSS)))
      } else { # Different endpoints
        # Note: this code will get not executed if nEps is 1.
        if(nRowTrt == nColTrt) { # Same treatment
          mCorrSim[p, q] <- mCorrSim[q, p] <- dEPCorr
        } else { # Different treatments
          mCorrSim[p, q] <- mCorrSim[q, p] <-
            sqrt((nRowSS/(nRowSS+nCtrlSS))*(nColSS/(nColSS+nCtrlSS))) * dEPCorr
        }
      }
    }
  }

  return(mCorrSim)
}

# Function for flattening the input of gMCP::calcPower() function into a dataframe
# row
FlattenGMCPInp <- function(inpCalPow) {
  # out <- calcPower(graph=gGraph, alpha=dAlpha, mean=dZMeans, corr.sim = mCorrSim,
  #                  corr.test = mCorrTest, n.sim = nSim)

  gGraph <- inpCalPow$graph
  mCorrSim <- inpCalPow$corr.sim
  mCorrTest <- inpCalPow$corr.test
  dZMeans <- inpCalPow$mean

  sG <- paste0("matrix(c(", paste0(as.vector(t(gGraph@m)), collapse = ","),
                  "), byrow=T, nrow=", nrow(gGraph@m), ")")

  sWt <- paste0("c(", paste0(gGraph@weights, collapse = ","), ")")

  sCorrSim <- paste0("matrix(c(", paste0(as.vector(t(mCorrSim)), collapse = ",")
                     , "), byrow=T, nrow=", nrow(mCorrSim), ")")

  sCorrTest <- if(is.null(mCorrTest)) {
    "NULL"
  } else {
    paste0("matrix(c(", paste0(as.vector(t(mCorrTest)), collapse = ",")
           , "), byrow=T, nrow=", nrow(mCorrTest), ")")
  }

  sMeans <- paste0("c(", paste0(dZMeans, collapse = ","), ")")

  dfOut <- cbind(sG, sWt, inpCalPow$alpha, sMeans, sCorrSim, sCorrTest,
                 inpCalPow$n.sim)

  colnames(dfOut) <- c("gMCP_G", "gMCP_Wts", "gMCP_alpha", "gMCP_means",
                       "gMCP_corr.sim", "gMCP_corr.test", "gMCP_n.sim")

  return(dfOut)
}

# A dummy function to replace calcPower() with zeroed output
DummyCalcPower <- function(graph, alpha, mean, corr.sim, corr.test, n.sim) {
  # Get the number of hypotheses from the graph matrix dimensions
  nHypotheses <- nrow(graph@m)

  # Create hypothesis names (H1, H2, etc.) for the LocalPower matrix
  hypothesisNames <- paste0("H", 1:nHypotheses)

  # Create a custom "out" object with all values set to 0
  localPowerMatrix <- matrix(0, nrow = nHypotheses, ncol = 1)
  rownames(localPowerMatrix) <- hypothesisNames

  out <- list(
    ExpRejections = 0,
    PowAtlst1 = 0,
    LocalPower = localPowerMatrix
  )

  return(out)
}

###################################
# test code
# nArms <- 5
# nTrtms <- nArms-1
# nEps <- 3
#
# for (p in 1:(nEps*nTrtms-1)) {
#   for (q in (p+1):(nEps*nTrtms)) {
#     nRowEP <- ceiling(p/nTrtms)
#     nRowTrt <- p - nTrtms*(nRowEP-1)
#
#     nColEP <- ceiling(q/nTrtms)
#     nColTrt <- q - nTrtms*(nColEP-1)
#
#     message(paste(p, nRowEP, nRowTrt, q, nColEP, nColTrt))
#   }
# }
#
