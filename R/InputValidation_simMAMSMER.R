# To test the input arguments of simMAMSMEP
# Note : The error text must contain the word "Invalid" so that later they can be processed appropriately
valInpsimMAMSMEP <- function(inps) {
  logs <- list()

  logs[[1]] <- ifelse(inps$Method == "CombPValue" || inps$Method == "CER", 0, "Invalid argument in 'Method'")

  logs[[2]] <- ifelse(inps$TestStatCont == "z" || inps$TestStatCont == "t-equal" || inps$TestStatCont == "t-unequal",
    0, "Invalid argument in 'TestStatCont'"
  )
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

  if(length(inps$InfoFrac) > 1){
    logs[[15]] <- ifelse(inps$SelectionScale == "delta" || inps$SelectionScale == "teststat" || inps$SelectionScale == "stderror" || inps$SelectionScale == "pvalue",
                         0, "Invalid argument in 'SelectionScale'"
    )

    logs[[16]] <- ifelse(inps$SelectionCriterion == "best" || inps$SelectionCriterion == "threshold" || inps$SelectionCriterion == "epsilon",
                         0, "Invalid argument in 'SelectionCriterion'"
    )

    logs[[18]] <- ifelse(inps$ImplicitSSR == "Selection" || inps$ImplicitSSR == "All" || inps$ImplicitSSR == "None",
                         0, "Invalid argument in 'ImplicitSSR'"
    )
  }else{
    # Fixed sample
    logs[[15]] <- 0
    logs[[16]] <- 0
    logs[[18]] <- 0
  }

  logs[[17]] <- ifelse(inps$Seed == "Random" || is.numeric(inps$Seed),
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


  logs[[21]] <- ifelse(inps$TestStatBin == "UnPooled" || inps$TestStatBin == "Pooled",
    0, "Invalid argument in 'TestStatBin'"
  )

  if(inps$nEps > 1){
    logs[[22]] <- ifelse(matrixcalc::is.positive.semi.definite(inps$EP.Corr),
                         0, "Invalid argument in 'EP.Corr', the matrix is not positive semi-definite")
  }



  return(logs)
}
