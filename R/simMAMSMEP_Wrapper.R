#' Wrapper function that takes all inputs through a single R dataframe
#' @param InputDF R Dataframe: This is the csv/excel input data in the R dataframe format
#' @export
simMAMSMEP_Wrapper <- function(InputDF){
  # Update the dataframe column names in the following mapping in case
  # the names in the input csv/excel changes
  lOut <- list()
  for (nModelNum in 1:nrow(InputDF)){
    # mapping to link simMAMSMEP function arguments with csv columns
    Method <- InputDF$Method[nModelNum]
    SampleSize <- InputDF$SampleSize[nModelNum]
    alpha <- InputDF$alpha[nModelNum]
    TestStat <- InputDF$TestStat[nModelNum]
    FWERControl <- InputDF$FWERControl[nModelNum]
    nArms <- InputDF$nArms[nModelNum]
    nEps <- InputDF$nEps[nModelNum]
    Arms.Mean <- eval(parse(text = InputDF$Arms.Mean[nModelNum]))
    Arms.std.dev <- eval(parse(text = InputDF$Arms.std.dev[nModelNum]))
    Arms.alloc.ratio <- eval(parse(text = InputDF$Arms.alloc.ratio[nModelNum]))
    EP.Corr <- eval(parse(text = InputDF$EP.Corr[nModelNum]))
    WI <- eval(parse(text = InputDF$WI[nModelNum]))
    G <- eval(parse(text = InputDF$G[nModelNum]))
    test.type <- InputDF$test.type[nModelNum]
    info_frac <- eval(parse(text = InputDF$info_frac[nModelNum]))
    typeOfDesign <- InputDF$typeOfDesign[nModelNum]
    MultipleWinners <- InputDF$MultipleWinners[nModelNum]
    Selection <- InputDF$Selection[nModelNum]
    SelectionLook <- InputDF$SelectionLook[nModelNum]
    SelectEndPoint <- InputDF$SelectEndPoint[nModelNum]
    SelectionScale <- InputDF$SelectionScale[nModelNum]
    SelectionCriterion <- InputDF$SelectionCriterion[nModelNum]
    SelectionParmeter <- InputDF$SelectionParmeter[nModelNum]
    KeepAssosiatedEps <- InputDF$KeepAssosiatedEps[nModelNum]
    ImplicitSSR <- InputDF$ImplicitSSR[nModelNum]
    nSimulation <- InputDF$nSimulation[nModelNum]
    Seed <- InputDF$Seed[nModelNum]
    SummaryStat <- InputDF$SummaryStat[nModelNum]
    plotGraphs <- InputDF$plotGraphs[nModelNum]
    Parallel <- InputDF$Parallel[nModelNum]
    # put the following code in try catch so the loop continues even if one iteration fails
    out <- simMAMSMEP(alpha = alpha, SampleSize = SampleSize,nArms = nArms,nEps = nEps, TestStat=TestStat, FWERControl = FWERControl,

                      Arms.Mean = Arms.Mean, Arms.std.dev = Arms.std.dev, Arms.alloc.ratio = Arms.alloc.ratio,

                      EP.Corr = EP.Corr,WI = WI, G = G, test.type = test.type,info_frac = info_frac,

                      typeOfDesign = typeOfDesign, MultipleWinners = MultipleWinners,

                      Selection = Selection,SelectionLook = SelectionLook,SelectEndPoint = SelectEndPoint,SelectionScale = SelectionScale,

                      SelectionCriterion = SelectionCriterion, SelectionParmeter = SelectionParmeter, KeepAssosiatedEps = KeepAssosiatedEps,

                      ImplicitSSR = ImplicitSSR, nSimulation = nSimulation, Seed = Seed, SummaryStat = SummaryStat,

                      Method = Method, plotGraphs = plotGraphs, Parallel = Parallel)
    # extract the power table from the output
    dfOverall_Powers_long <- out$Overall_Powers_df
    # convert the power table to wide format and add serial number column
    dfOverall_Powers_wide <- dfOverall_Powers_long %>%
      pivot_wider(names_from = Overall_Powers, values_from = Values) %>%
      mutate(Sno = nModelNum) %>%
      relocate(Sno, .before = everything())

    # add each iterations power table to a list
    lOut[[nModelNum]] <- dfOverall_Powers_wide
    print(paste0("Power table for ", nModelNum, ":"))
    print(dfOverall_Powers_wide)
    print(paste0("Model ", nModelNum, " execution completed successfully."))
  }
  # rbind power tables for each iteration to produce a single table
  dfOut <- do.call(rbind, lOut)
  return(dfOut)
}
