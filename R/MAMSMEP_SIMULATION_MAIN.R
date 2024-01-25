#' Function to perform Adaptive GMCP simulation for Multi-Arm Multi-Stage Multi-Endpoint simulations for Combining p-values method and CER method(2-Stage)
#' @param alpha Type-1 error
#' @param SampleSize integer valued Sample Size(default: 500)
#' @param nArms integer value to specify the number of arms (default: 3)
#' @param nEps integer value to specify the number of endpoints
#' @param Arms.Mean Numeric list to specify the arm-wise mean for each endpoint
#' @param Arms.std.dev Numeric list to specify the arm-wise standard deviation for each endpoint
#' @param Arms.alloc.ratio Numeric Vector to specify the arm-wise allocation ratio
#' @param EP.Corr correlation matrix for the endpoints(Normal)
#' @param WI Vector of Initial Weights for Global Null.
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
#' @example ./internalData/MAMSMEP_Simulation_Example.R
#' @export
simMAMSMEP <- function(
    alpha = 0.025,
    SampleSize = 500,
    nArms = 3,
    nEps  = 2,
    Arms.Mean = list('EP1' = c(0,0.4,0.3),
                     'EP2' = c(0.1,0.45,0.25)),
    Arms.std.dev = list('EP1' = c(1,1,1), 'EP2' = c(1,1,1)),
    Arms.alloc.ratio = c(1,1,1),
    EP.Corr = matrix(c(1,0.5,
                       0.5,1),
                     nrow = 2),
    WI =  c(0.5,0.5,0,0),
    G = matrix(c(0,0.5,0.5,0,
                 0.5,0,0,0.5,
                 0,1,0,0,
                 1,0,0,0),
               nrow = nEps*(nArms-1), byrow = T),
    test.type = 'Partly-Parametric',
    info_frac = c(0.5,1),
    typeOfDesign = "asOF",
    MultipleWinners = TRUE,
    Selection= TRUE,
    SelectionLook = 1,
    SelectEndPoint = 1,
    SelectionScale = 'pvalue',
    SelectionCriterion ='best',
    SelectionParmeter = 1,
    KeepAssosiatedEps = TRUE,
    ImplicitSSR = 'All',
    nSimulation = 1000,
    Seed = 100,
    SummaryStat = FALSE,
    Method = 'CombPValue',
    plotGraphs = TRUE
)
{
  TailType <- 'RightTail'       ##Default Right
  UpdateStrategy <- F           ##Not implemented yet
  des.type <- 'MAMSMEP'         ##Multi-Arm Multi-Stage Multi-EndPoints

  #Object to run Simulations
  gmcpSimObj <- list(
    #Methodology
    'Method' = Method,

    #test parameters
    'nArms'= nArms,                            'nLooks'= length(info_frac),
    'nEps'= nEps,                              'nHypothesis' = nEps*(nArms-1),
    'TailType' = TailType,                     'des.type' = des.type,
    'Max_SS' = SampleSize,                     'test.type' = test.type,
    'IntialWeights'=WI,                        'G' = G,
    'Correlation' = NA,                        'alpha' = alpha,

    #Boundary
    'InfoFrac' = info_frac,                    'typeOfDesign'=typeOfDesign,

    #Multiple Winners
    'MultipleWinners'=MultipleWinners,

    #Response Generation
    'Arms.Mean' = Arms.Mean,                   'Arms.std.dev' = Arms.std.dev,
    'Arms.alloc.ratio' = Arms.alloc.ratio,     'Arms.alloc.ratio' = Arms.alloc.ratio,
    'EP.Corr' = EP.Corr,

    #Selection
    'SelectEndPoint'= SelectEndPoint,         'Selection' = Selection,
    'SelectionLook' = SelectionLook,          'SelectionScale' = SelectionScale,
    'SelectionCriterion' = SelectionCriterion,'SelectionParmeter' = SelectionParmeter,
    'KeepAssosiatedEps' = KeepAssosiatedEps,

    #Simulation Parameters
    'nSimulation' = nSimulation,               'Seed' = Seed,
    'SummaryStat' = SummaryStat,

    #SSR
    'ImplicitSSR' = ImplicitSSR,

    #Update Strategy
    'UpdateStrategy' =UpdateStrategy,

    #Graph Plot
    'plotGraphs' = plotGraphs
  )

  out <- mnMAMSMEP_sim2(gmcpSimObj)
  if(plotGraphs){
    #out$iniGraph
    return(out$DetailOutTabs)
  }else
  {
    return(out$DetailOutTabs)
  }
}
