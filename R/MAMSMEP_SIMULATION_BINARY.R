# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

simMAMSMEP_BIN <- function(
    Method = "CombPValue",
    alpha = 0.025,
    SampleSize = 500,
    TestStat = "t-equal",
    FWERControl = "CombinationTest",
    nArms = 3,
    nEps = 2,
    Arms.Mean = list(
      "EP1" = c(0, 0.4, 0.3),
      "EP2" = c(0.1, 0.45, 0.25)
    ),
    Arms.std.dev = list("EP1" = c(1, 1, 1), "EP2" = c(1, 1, 1)),
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
    "TestStat" = TestStat, "alpha" = alpha,
    "nArms" = nArms, "nLooks" = length(info_frac),
    "nEps" = nEps, "nHypothesis" = nEps * (nArms - 1),
    "TailType" = TailType, "des.type" = des.type,
    "Max_SS" = SampleSize, "test.type" = test.type,
    "IntialWeights" = WI, "G" = G,
    "Correlation" = NA, "FWERControl" = FWERControl,

    # Boundary
    "InfoFrac" = info_frac, "typeOfDesign" = typeOfDesign,

    # Multiple Winners
    "MultipleWinners" = MultipleWinners,

    # Response Generation
    "Arms.Mean" = Arms.Mean, "Arms.std.dev" = Arms.std.dev,
    "Arms.alloc.ratio" = Arms.alloc.ratio, "Arms.alloc.ratio" = Arms.alloc.ratio,
    "EP.Corr" = EP.Corr,

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
  out <- mnMAMSMEP_sim2(gmcpSimObj)
  if (plotGraphs) {
    # out$iniGraph
    return(out$DetailOutTabs)
  } else {
    return(out$DetailOutTabs)
  }
}
