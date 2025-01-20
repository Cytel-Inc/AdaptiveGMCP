# The file contains supporting functions(Detailed Output Tables) for Adaptive GMCP Simulation Function AdaptGmcp_Simulation/adaptGMCP_SIM(.) ----
## Author: Ajoy.M

#-------------- -
# Compute Overall Power Table from all the simulation
#-------------- -
SimPowers <- function(nSimulation, nSimulation_Stage2, PowerTab) {
  values <- as.numeric(apply(PowerTab[, -1], 2, function(x) {
    sum(x) / (nrow(PowerTab))
  }))
  # 95% confidence interval
  z_alpha <- qnorm(1-0.025)
  UL <- values + z_alpha*sqrt(values*(1-values)/(nrow(PowerTab)))
  LL <- values - z_alpha*sqrt(values*(1-values)/(nrow(PowerTab)))
  ConfIntv <- sapply(1:length(values),function(i){
    paste('(',round(LL[i],5),',',round(UL[i],5),')',sep = '')
  })

  PowerTable <- data.frame(
    "Overall_Powers" = c("Global Power", "Conjunctive Power", "Disjunctive Power", "FWER"),
    "Values" = values,
    "ConfIntv_95perc" = ConfIntv
  )
  PowerTable
}

#------------- -
# Count contribution to different powers from each simulations
#------------- -
CountPower <- function(simID, SummaryStatFile, TrueNull) {
  rejMat <- subset(SummaryStatFile, SimID == simID)
  rejMat <- rejMat[, grep("RejStatus", names(rejMat))]
  rej.final <- apply(rejMat, 2, function(col) any(col, na.rm = T))
  data.frame(
    "simID" = simID,
    # "simID_Stage2" = simID_Stage2,
    "nG" = as.integer(any(rej.final)),
    "nC" = as.integer(
      ifelse(length(rej.final[!TrueNull]) == 0, 0, all(rej.final[!TrueNull]))
    ),
    "nD" = as.integer(
      ifelse(length(rej.final[!TrueNull]) == 0, 0, any(rej.final[!TrueNull]))
    ),
    "nF" = as.integer(
      ifelse(length(rej.final[TrueNull]) == 0, 0, any(rej.final[TrueNull]))
    )
  )
}

#------------ -
# Identify True Null based of the response generation inputs
#------------ -
checkTrueNull <- function(gmcpSimObj, index_map) {
  pairsHyp <- strsplit(index_map$Names, split = "-")

  delta <- unlist(lapply(pairsHyp, function(x) {
    abs(gmcpSimObj$Arms.Mean[x[1]] - gmcpSimObj$Arms.Mean[x[2]])
  }))
  index_map$TrueNull <- (delta < 1E-6)
  index_map
}

#------------ -
# Identify True Null based of the response generation inputs
#------------ -
checkTrueNull2 <- function(HypoMap, Arms.Mean) {
  delta <- unlist(lapply(1:nrow(HypoMap), function(i) {
    abs(Arms.Mean[[HypoMap$Groups[i]]][HypoMap$Treatment[i]] -
      Arms.Mean[[HypoMap$Groups[i]]][HypoMap$Control[i]])
  }))
  TrueNull <- (delta < 1E-6)
  TrueNull
}

#------------ -
# Identify True Null based of the response generation inputs
#------------ -
checkTrueNull3 <- function(HypoMap, Arms.Mean, Arms.Prop) {
  delta <- unlist(lapply(1:nrow(HypoMap), function(i) {
    if (HypoMap$EpType[i] == "Continuous") {
      abs(Arms.Mean[[HypoMap$Groups[i]]][HypoMap$Treatment[i]] -
        Arms.Mean[[HypoMap$Groups[i]]][HypoMap$Control[i]])
    } else if (HypoMap$EpType[i] == "Binary") {
      abs(Arms.Prop[[HypoMap$Groups[i]]][HypoMap$Treatment[i]] -
        Arms.Prop[[HypoMap$Groups[i]]][HypoMap$Control[i]])
    }
  }))
  TrueNull <- (delta < 1E-6)
  TrueNull
}

#------------- -
# Count contribution to different powers from each simulations
#------------- -
CountEfficacy <- function(simID, SummaryStatFile, interHypo = NULL) {
  # Subset rows and select relevant columns using base R for better parallel compatibility
  rejMat <- SummaryStatFile[SummaryStatFile$SimID == simID,
                            grep("RejStatus", colnames(SummaryStatFile), value = TRUE)]

  # Check rejection status for each column
  rej.final <- colSums(!is.na(rejMat) & rejMat != 0) > 0
  n <- length(rej.final)
  if (n == 0) rej.final <- rep(FALSE, n)

  # if (is.null(interHypo)){
  # Generate combinations - moved this calculation to PreSimObj for CER
  interHypo <- genCombs(n)
  # }

  # Create index strings with explicit type conversion
  idx <- apply(interHypo, 1, function(x) paste(as.character(x), collapse = ""))
  rej_idx <- paste(as.character(as.integer(rej.final)), collapse = "")

  # Generate column names
  col_name <- apply(interHypo, 1, function(row) {
    paste(paste0("H", which(row == 1)), collapse = ",")
  })

  # Mark efficacy
  eff <- numeric(length(idx))
  match_idx <- match(rej_idx, idx)
  if (!is.na(match_idx)) {
    eff[match_idx] <- 1
  }

  # Create output data frame using base R
  eff_count <- data.frame(
    simID = simID,
    matrix(eff, nrow = 1)
  )
  colnames(eff_count) <- c("simID",col_name)

  return(eff_count)
}

genCombs <- function(n) {
  if (n == 0) {
    return(data.frame())
  }
  # Calculate total number of combinations (2^n - 1, excluding all zeros)
  total_combs <- 2^n - 1

  # Create binary representations
  binary_nums <- 1:total_combs

  # Convert to binary matrix
  result <- matrix(0, nrow = total_combs, ncol = n)
  for(i in 1:n) {
    result[, i] <- (binary_nums %% 2^i) >= 2^(i-1)
  }

  result <- as.data.frame(result)
  rownames(result) <- NULL
  return(result)
}
