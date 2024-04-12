#' Function to perform Adaptive GMCP simulation for Multi-Arm Multi-Stage Multi-Endpoint simulations for Combining p-values method and CER method(2-Stage)
#' @param Method 'CombPValue': for combining p-values method, 'CER': for Conditional Error method.
#' @param alpha Type-1 error
#' @param SampleSize integer valued Sample Size(default: 500)
#' @param TestStatCont Test Statistics for continuous endpoints; options: 't-equal' : for t statistics with equal variance, 't-unequal' : for t statistics with unequal variance, 'z' for z statistics
#' @param TestStatBin Test Statistics for Binary endpoints; options: 'UnPooled', 'Pooled'
#' @param FWERControl applicable for CER method only, 'CombinationTest': combined two stage incremental test statistics, 'None': Cumulative test statistics.
#' @param nArms integer value to specify the number of arms (default: 3)
#' @param nEps integer value to specify the number of endpoints
#' @param lEpType list with endpoint types
#' @param Arms.Mean Numeric list to specify the arm-wise mean for each endpoint; Note: The first input is for control arm and the rest are for the treatments.
#' @param Arms.std.dev Numeric list to specify the arm-wise standard deviation for each endpoint; Note: The first input is for control arm and the rest are for the treatments.
#' @param CommonStdDev TRUE = the treatment standard deviations assumed to be same as the control for boundary computations for continuous endpoints, FALSE = the treatment standard deviations assumed to be same as given in Arms.std.dev.
#' @param Arms.Prop Numeric list to specify the arm-wise proportions for each endpoint; Note: The first input is for control arm and the rest are for the treatments.
#' @param Arms.alloc.ratio Numeric Vector to specify the arm-wise allocation ratio; Note: The first input is for control arm and the rest are for the treatments.
#' @param EP.Corr correlation matrix for the endpoints(Normal)
#' @param WI Vector of Initial Weights for Global Null; Note: Hypotheses will follow the order of Endpoints and Treatments as given in 'Arms.Mean' and 'Arms.std.dev' inputs e.g.: If 'Arms.Mean' are given in the format list('EP1'=c(ctr_mean, trt1_mean, trt2_mean), 'EP2'=c(ctr_mean, trt1_mean, trt2_mean)) then the four hypotheses will be H1 = (Trt1 vs Ctr for EP1), H2 = (Trt2 vs Ctr for EP1), H3 = (Trt1 vs Ctr for EP2), H4 = (Trt2 vs Ctr for EP2), The initial weights and the transition matrix will follow the order of hypothesis accordingly as (H1,H2,H3,H4)
#' @param G  Numeric Matrix to specify the Transition Matrix.
#' @param test.type Character to specify the type of test want to perform; Available tests for Combining P-values Method :- 'Bonf': Bonferroni, 'Sidak': Sidak, 'Simes': Simes, 'Dunnett': Dunnett and  'Partly-Parametric': Mixed type Tests. Available tests for CER Method :- "Parametric": Weighted Dunnett , "Non-Parametric": Weighted Bonferroni and  'Partly-Parametric': Mixed type Tests.
#' @param info_frac Numeric Vector to specify look position as fraction of sample size.(for one look can be specified as 1)
#' @param typeOfDesign The type of design. Type of design is one of the following: O'Brien & Fleming ("OF"), Pocock ("P"), Wang & Tsiatis Delta class ("WT"), Pampallona & Tsiatis ("PT"), Haybittle & Peto ("HP"), Optimum design within Wang & Tsiatis class ("WToptimum"), O'Brien & Fleming type alpha spending ("asOF"), Pocock type alpha spending ("asP"), Kim & DeMets alpha spending ("asKD"), Hwang, Shi & DeCani alpha spending ("asHSD"), no early efficacy stop ("noEarlyEfficacy"), default is "OF".
#' @param MultipleWinners Logical; TRUE: Stop the trial only no more efficacy is possible, FALSE: Stop if at-least one efficacy is observed
#' @param Selection Logical: TRUE if selection required at interim(default = FALSE)
#' @param SelectionLook Numeric Vector to specify the selection looks
#' @param SelectEndPoint Indicator to specify which endpoint to select from, e.g. '1': Endpoint 1, '2':Endpoint 2, 'overall': overall
#' @param SelectionScale Character: Scale parameter on which selection will be based on, options 'delta': delta, 'teststat': Test Statistics, 'stderror' : Standard Error of the test stat,  'pvalue': p-value(un-adj) based selection
#' @param SelectionCriterion Character: 'best': best r, 'threshold': threshold for selection, 'epsilon': for epsilon neighborhood
#' @param SelectionParmeter r for best, threshold value for threshold or epsilon distance
#' @param KeepAssosiatedEps Logical, True: keep all the associated hypothesis for the selected arms
#' @param ImplicitSSR Character; 'Selection': re-allocate samples only from de-selected arms to available arms, 'All': Allocate all the planned samples(for the look) to the available arms, 'None': No Re-allocation
#' @param UpdateStrategy Logical to specify the updated strategy (Not Implemented yet) default FALSE
#' @param nSimulation Numeric: number of simulations(default=1000)
#' @param Seed 'Random' for randomly generating seed else any integer value(default = 'Random')
#' @param SummaryStat Logical; TRUE if simulation level data is required(default = FALSE)
#' @param plotGraphs Logical; TRUE: plot the initial graph
#' @param Parallel Logical; TRUE: Parallel computations
#' @example ./internalData/MAMSMEP_Simulation_Example.R
#' @export
simMAMSMEP <- function(
    Method = "CombPValue",
    alpha = 0.025,
    SampleSize = 500,
    TestStatCont = "t-equal",
    TestStatBin = "UnPooled",
    FWERControl = "None",
    nArms = 3,
    nEps = 2,
    lEpType = list(
      "EP1" = "Continuous",
      "EP2" = "Binary"
    ),
    Arms.Mean = list(
      "EP1" = c(0, 0.4, 0.3),
      "EP2" = NA
    ),
    Arms.std.dev = list(
      "EP1" = c(1.1, 1.2, 1.3),
      "EP2" = NA
    ),
    CommonStdDev = FALSE,
    Arms.Prop = list(
      "EP1" = NA,
      "EP2" = c(0.2, 0.35, 0.45)
    ),
    Arms.alloc.ratio = c(1, 1, 1),
    EP.Corr = matrix(
      c(
        1, 0.5,
        0.5, 1
      ),
      nrow = 2
    ),
    WI = c(0.5, 0.5, 0, 0),
    G = matrix(
      c(
        0, 0.5, 0.5, 0,
        0.5, 0, 0, 0.5,
        0, 1, 0, 0,
        1, 0, 0, 0
      ),
      nrow = nEps * (nArms - 1), byrow = T
    ),
    test.type = "Partly-Parametric",
    info_frac = c(0.5, 1),
    typeOfDesign = "asOF",
    MultipleWinners = TRUE,
    Selection = TRUE,
    SelectionLook = 1,
    SelectEndPoint = 1,
    SelectionScale = "pvalue",
    SelectionCriterion = "best",
    SelectionParmeter = 1,
    KeepAssosiatedEps = TRUE,
    ImplicitSSR = "All",
    nSimulation = 1000,
    Seed = 100,
    SummaryStat = FALSE,
    plotGraphs = TRUE,
    Parallel = TRUE) {
  Parallel <- Parallel
  TailType <- "RightTail" ## Default Right
  UpdateStrategy <- F ## Not implemented yet
  des.type <- "MAMSMEP" ## Multi-Arm Multi-Stage Multi-EndPoints


  # Object to run Simulations
  gmcpSimObj <<- list(
    # Methodology
    "Method" = Method,

    # test parameters
    "TestStatCont" = TestStatCont, "alpha" = alpha,
    "nArms" = nArms, "nLooks" = length(info_frac),
    "nEps" = nEps, "nHypothesis" = nEps * (nArms - 1),
    "TailType" = TailType, "des.type" = des.type,
    "Max_SS" = SampleSize, "test.type" = test.type,
    "IntialWeights" = WI, "G" = G,
    "Correlation" = NA, "FWERControl" = FWERControl,
    "lEpType" = lEpType, "TestStatBin" = TestStatBin,

    # Boundary
    "InfoFrac" = info_frac, "typeOfDesign" = typeOfDesign,
    "CommonStdDev" = CommonStdDev,

    # Multiple Winners
    "MultipleWinners" = MultipleWinners,

    # Response Generation
    "Arms.Mean" = Arms.Mean, "Arms.std.dev" = Arms.std.dev,
    "Arms.alloc.ratio" = Arms.alloc.ratio, "Arms.alloc.ratio" = Arms.alloc.ratio,
    "EP.Corr" = EP.Corr, "Arms.Prop" = Arms.Prop,
    "prop.ctr" = lapply(Arms.Prop, function(x) {
      x[2]
    }),

    # Selection
    "SelectEndPoint" = SelectEndPoint, "Selection" = Selection,
    "SelectionLook" = SelectionLook, "SelectionScale" = SelectionScale,
    "SelectionCriterion" = SelectionCriterion, "SelectionParmeter" = SelectionParmeter,
    "KeepAssosiatedEps" = KeepAssosiatedEps,

    # Simulation Parameters
    "nSimulation" = nSimulation, "Seed" = Seed,
    "SummaryStat" = SummaryStat, "Parallel" = Parallel,

    # SSR
    "ImplicitSSR" = ImplicitSSR,

    # Update Strategy
    "UpdateStrategy" = UpdateStrategy,

    # Graph Plot
    "plotGraphs" = plotGraphs
  )

  logs <- valInpsimMAMSMEP(inps = gmcpSimObj)
  FailedLogs <- logs[sapply(logs, function(lg) lg != 0)]

  if (length(FailedLogs) != 0) {
    return(FailedLogs)
    stop("Input Error")
  }
  out <- MAMSMEP_sim2(gmcpSimObj)
  if (plotGraphs) {
    # out$iniGraph
    return(out$DetailOutTabs)
  } else {
    return(out$DetailOutTabs)
  }
}
