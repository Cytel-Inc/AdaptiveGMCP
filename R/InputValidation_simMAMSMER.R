# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

# To test the input arguments of simMAMSMEP
# Note : The error text must contain the word "Invalid" so that later they can be processed appropriately
valInpsimMAMSMEP <- function(inps) {
  logs <- list()

  logs[[1]] <- ifelse(inps$Method == "CombPValue" || inps$Method == "CER", 0, "Invalid argument in 'Method'")

  if(all(inps$lEpType == 'Continuous')){
    logs[[2]] <- ifelse(inps$TestStatCont == "z" || inps$TestStatCont == "t-equal" || inps$TestStatCont == "t-unequal",
                        0, "Invalid argument in 'TestStatCont'"
    )
  }else{
    logs[[2]] <- 0
  }

  logs[[3]] <- ifelse(inps$nArms >= 2, 0, "Invalid argument in 'nArms'")

  logs[[4]] <- ifelse(ifelse(inps$Method == "CER", length(inps$InfoFrac) <= 2 & length(inps$InfoFrac) >= 1, length(inps$InfoFrac) >= 1),
    0, "Invalid argument in 'info_frac'"
  )

  logs[[5]] <- ifelse(inps$nEps >= 1, 0, "Invalid argument in 'nEps'")

  logs[[6]] <- ifelse(inps$Max_SS > 0 & is.numeric(inps$Max_SS), 0, "Invalid argument in 'SampleSize'")

  logs[[7]] <- ifelse(ifelse(inps$Method == "CER", inps$test.type == "Parametric" || inps$test.type == "Partly-Parametric" || inps$test.type == "Non-Parametric",
    inps$test.type == "Partly-Parametric" || inps$test.type == "Dunnett" || inps$test.type == "Bonf" || inps$test.type == "Sidak" || inps$test.type == "Simes"
  ),
  0, "Invalid argument in 'test.type'"
  )

  logs[[8]] <- ifelse(length(inps$IntialWeights) == inps$nHypothesis & sum(inps$IntialWeights) <= 1,
    0, "Invalid argument in 'WI'"
  )

  logs[[9]] <- ifelse(all(diag(inps$G) == 0) & all(rowSums(inps$G) <= 1) & nrow(inps$G) == inps$nHypothesis & ncol(inps$G) == inps$nHypothesis,
    0, "Invalid argument in 'G'"
  )

  logs[[10]] <- ifelse(ifelse(inps$Method == "CER", inps$FWERControl == "CombinationTest" || inps$FWERControl == "None", TRUE),
    0, "Invalid argument in 'FWERControl'"
  )

  if(all(inps$lEpType == 'Continuous')){
    logs[[11]] <- ifelse(length(inps$Arms.Mean) == inps$nEps,
                         0, "Invalid argument in 'Arms.Mean'"
    )
    logs[[12]] <- ifelse(length(inps$Arms.std.dev) == inps$nEps,
                         0, "Invalid argument in 'Arms.std.dev'"
    )
  }else{
    logs[[11]] <- 0
    logs[[12]] <- 0
  }


  logs[[13]] <- ifelse(length(inps$Arms.alloc.ratio) == inps$nArms,
    0, "Invalid argument in 'Arms.alloc.ratio'"
  )

  logs[[14]] <- ifelse(ifelse(inps$nEps >= 2, nrow(inps$EP.Corr) == inps$nEps & ncol(inps$EP.Corr) == inps$nEps, TRUE),
    0, "Invalid argument in 'EP.Corr'"
  )

  if(length(inps$InfoFrac) > 1 & inps$Selection){
    logs[[15]] <- ifelse(inps$SelectionScale == "delta" || inps$SelectionScale == "teststat" ||
                           inps$SelectionScale == "stderror" || inps$SelectionScale == "pvalue" || inps$SelectionCriterion == "random",
                         0, "Invalid argument in 'SelectionScale'"
    )

    logs[[16]] <- ifelse(inps$SelectionCriterion == "best" || inps$SelectionCriterion == "threshold" ||
                           inps$SelectionCriterion == "epsilon" || inps$SelectionCriterion == "random",
                         0, "Invalid argument in 'SelectionCriterion'"
    )
  }else{
    # Fixed sample
    logs[[15]] <- 0
    logs[[16]] <- 0
  }
  if(length(inps$InfoFrac) > 1){
    logs[[18]] <- ifelse(inps$ImplicitSSR == "Selection" || inps$ImplicitSSR == "All" || inps$ImplicitSSR == "None",
                         0, "Invalid argument in 'ImplicitSSR'"
    )
  }else{
    logs[[18]] <- 0
  }

  logs[[17]] <- ifelse(inps$Seed == "Random" || is.numeric(as.numeric(inps$Seed)),
    0, "Invalid argument in 'Seed'"
  )



  logs[[19]] <- ifelse(length(inps$lEpType) == inps$nEps,
    0, "Invalid argument in 'lEpType'"
  )

  if(all(inps$lEpType == 'Binary')){
    logs[[20]] <- ifelse(length(inps$Arms.Prop) == inps$nEps,
                         0, "Invalid argument in 'Arms.Prop'"
    )
  }else{
    logs[[20]] <- 0
  }

  if(all(inps$lEpType == 'Binary')){
    logs[[21]] <- ifelse(inps$TestStatBin == "UnPooled" || inps$TestStatBin == "Pooled",
                         0, "Invalid argument in 'TestStatBin'"
    )
  }else{
    logs[[21]] <- 0
  }


  if(inps$nEps > 1){
    logs[[22]] <- ifelse(matrixcalc::is.positive.semi.definite(inps$EP.Corr),
                         0, "Invalid argument in 'EP.Corr', the matrix is not positive semi-definite")
  }else {
     logs[[22]] <- 0
  }

  # ---- Survival endpoint validation ----
  isSurvival <- any(sapply(inps$lEpType, function(x) x == "Survival"))

  if (isSurvival) {
    # Survival requires exactly 1 endpoint
    logs[[23]] <- ifelse(inps$nEps == 1,
      0, "Invalid argument in 'nEps': must be 1 for Survival endpoint")

    # Only fixed sample and 2-look (1 interim + 1 final) designs supported for Survival
    logs[[24]] <- ifelse(inps$nLooks <= 2,
      0, "Invalid argument in 'info_frac': Only 1 or 2 looks are supported for Survival endpoint")

    # Only CombPValue method supported for Survival
    logs[[25]] <- ifelse(inps$Method == "CombPValue",
      0, "Invalid argument in 'Method': only 'CombPValue' is supported for Survival endpoint")

    # totalEvents must be specified, positive integer, and <= SampleSize
    if (!is.na(inps$totalEvents)) {
      logs[[26]] <- ifelse(
        is.numeric(inps$totalEvents) && length(inps$totalEvents) == 1 &&
          inps$totalEvents > 0 && inps$totalEvents == round(inps$totalEvents) &&
          inps$totalEvents <= inps$Max_SS,
        0, "Invalid argument in 'totalEvents': must be a positive integer <= SampleSize")
    } else {
      logs[[26]] <- "Invalid argument in 'totalEvents': must be specified for Survival endpoint"
    }

    # armsHazardRates must be specified for Survival endpoint(s).
    # Expected structure:
    # - list of length nEps
    # - for non-Survival endpoints: element must be NA (length 1)
    # - for Survival endpoints: element must be a positive numeric vector of length nArms
    surv.idx <- which(inps$lEpType == "Survival")
    ahr <- inps$armsHazardRates

    armsHazardRatesValid <- TRUE
    armsHazardRatesMsg <- "Invalid argument in 'armsHazardRates': must be a list of length nEps with NA for non-Survival endpoints and positive numeric vectors of length nArms for Survival endpoints"

    if (is.null(ahr)) {
      armsHazardRatesValid <- FALSE
      armsHazardRatesMsg <- "Invalid argument in 'armsHazardRates': must be specified for Survival endpoint(s)"
    } else if (!is.list(ahr)) {
      armsHazardRatesValid <- FALSE
    } else if (length(ahr) != inps$nEps) {
      armsHazardRatesValid <- FALSE
      armsHazardRatesMsg <- "Invalid argument in 'armsHazardRates': must be a list of length nEps"
    } else {
      for (i in seq_len(inps$nEps)) {
        elem <- ahr[[i]]
        if (i %in% surv.idx) {
          if (!(is.numeric(elem) &&
            length(elem) == inps$nArms &&
            !any(is.na(elem)) &&
            all(elem > 0))) {
            armsHazardRatesValid <- FALSE
            armsHazardRatesMsg <- paste0(
              "Invalid argument in 'armsHazardRates': element ", i,
              " (Survival endpoint) must be a positive numeric vector of length nArms"
            )
            break
          }
        } else {
          is.na.scalar <- length(elem) == 1 && is.atomic(elem) && is.na(elem)
          if (!is.na.scalar) {
            armsHazardRatesValid <- FALSE
            armsHazardRatesMsg <- paste0(
              "Invalid argument in 'armsHazardRates': element ", i,
              " (non-Survival endpoint) must be NA"
            )
            break
          }
        }
      }
    }

    logs[[27]] <- ifelse(armsHazardRatesValid, 0, armsHazardRatesMsg)

    # accrualStartTimes and accrualRates must be specified and same length
    accrualStartTimesValid <- !all(is.na(inps$accrualStartTimes)) &&
      is.numeric(inps$accrualStartTimes)
    accrualRatesValid <- !all(is.na(inps$accrualRates)) &&
      is.numeric(inps$accrualRates)

    if (accrualStartTimesValid && accrualRatesValid) {
      logs[[28]] <- ifelse(
        length(inps$accrualStartTimes) == length(inps$accrualRates),
        0, "Invalid argument in 'accrualStartTimes'/'accrualRates': must be same length")

      logs[[29]] <- ifelse(inps$accrualStartTimes[1] == 0,
        0, "Invalid argument in 'accrualStartTimes': first element must be 0")

      logs[[30]] <- ifelse(all(inps$accrualRates > 0),
        0, "Invalid argument in 'accrualRates': all values must be positive")

      # accrualStartTimes must be in strictly increasing order
      if (length(inps$accrualStartTimes) > 1) {
        logs[[31]] <- ifelse(
          all(diff(inps$accrualStartTimes) > 0),
          0, "Invalid argument in 'accrualStartTimes': must be in strictly increasing order")
      } else {
        logs[[31]] <- 0
      }
    } else {
      logs[[28]] <- "Invalid argument in 'accrualStartTimes'/'accrualRates': both must be specified for Survival endpoint"
      logs[[29]] <- 0
      logs[[30]] <- 0
      logs[[31]] <- 0
    }
  } else {
    # Not a survival trial: skip survival-specific checks
    logs[[23]] <- 0
    logs[[24]] <- 0
    logs[[25]] <- 0
    logs[[26]] <- 0
    logs[[27]] <- 0
    logs[[28]] <- 0
    logs[[29]] <- 0
    logs[[30]] <- 0
    logs[[31]] <- 0
  }

  return(logs)
}
