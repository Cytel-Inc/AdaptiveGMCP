#' #' Function to perform Adaptive GMCP simulation for Multi-Arm designs
#' #' @param alpha Type-1 error
#' #' @param SampleSize integer valued Sample Size(default: 200)
#' #' @param nArms integer value to specify the number of arms (default: 4)
#' #' @param Arms.Name Optional, character vector to contain the names of the arms(\code{e.g.- c('Placebo','ActiveCTR','D1','D2','D3','D4')}), defaults = NA
#' #' @param Arms.Mean Numeric Vector to specify the arm-wise mean for response generation
#' #' @param Arms.std.dev Numeric Vector to specify the arm-wise standard deviation for response generation
#' #' @param Arms.alloc.ratio Numeric Vector to specify the arm-wise allocation ratio
#' #' @param WI Vector of Initial Weights for Global Null.
#' #' @param G  Numeric Matrix to specify the Transition Matrix.
#' #' @param test.type Character to specify the type of test want to perform; 'Bonf': Bonferroni, 'Sidak': Sidak, 'Simes': Simes, 'Dunnett': Dunnett and  'Partly-Parametric': Partly Parametric Tests
#' #' @param Correlation Numeric Matrix to specify the correlation between test statistics, NA if the correlation is unknown. (applicable only for Dunnett and Partly-Parametric tests)
#' #' @param info_frac Numeric Vector to specify look position as fraction of sample size.(for one look can be specified as 1)
#' #' @param typeOfDesign The type of design. Type of design is one of the following: O'Brien & Fleming ("OF"), Pocock ("P"), Wang & Tsiatis Delta class ("WT"), Pampallona & Tsiatis ("PT"), Haybittle & Peto ("HP"), Optimum design within Wang & Tsiatis class ("WToptimum"), O'Brien & Fleming type alpha spending ("asOF"), Pocock type alpha spending ("asP"), Kim & DeMets alpha spending ("asKD"), Hwang, Shi & DeCani alpha spending ("asHSD"), no early efficacy stop ("noEarlyEfficacy"), default is "OF".
#' #' @param MultipleWinners Logical; TRUE: Stop the trial only no more efficacy is possible, FALSE: Stop if at-least one efficacy is observed
#' #' @param Selection Logical: TRUE if selection required at interim(default = FALSE)
#' #' @param SelectionLook Numeric Vector to specify the selection looks
#' #' @param SelectionMethods Character: Scale parameter on which selection will be based on; 'pvalue': p-value(un-adj) based selection, 'ArmID': Arm-Wise selection
#' #' @param SelectionCriterion Character: 'best': best r, 'threshold': threshold for selection
#' #' @param SelectionParmeter r for best, threshold value for threshold, Arms Index or Names for Arm-Wise selection
#' #' @param ImplicitSSR Character; 'Selection': re-allocate samples only from de-selected arms to available arms, 'All': Allocate all the planned samples(for the look) to the available arms, 'None': No Re-allocation
#' #' @param UpdateStrategy Logical to specify the updated strategy (Not Implemented yet) default FALSE
#' #' @param nSimulation Numeric: number of simulations(default=1000)
#' #' @param Seed 'Random' for randomly generating seed else any integer value(default = 'Random')
#' #' @param SummaryStat Logical; TRUE if simulation level data is required(default = FALSE)
#' adaptGMCP_SIM <- function(
#'     alpha = 0.025,
#'     SampleSize = 200,
#'     nArms = 4,
#'     Arms.Name = NA,
#'     Arms.Mean = c(0,0,0.1,0.4),
#'     Arms.std.dev = c(1,1,1,1),
#'     Arms.alloc.ratio = c(1,1,1,1),
#'     WI = c(0.2,0.3,0.5),
#'     G = matrix(c(0,0.25,0.75,
#'                  0.25,0,0.75,
#'                  0.25,0.75,0), byrow = T, nrow = 3)
#'     ,
#'     test.type = 'Partly-Parametric',
#'     Correlation = matrix(c(1,NA,NA,
#'                            NA,1,0.5,
#'                            NA,0.5,1), byrow = T, nrow = 3),
#'     info_frac = c(0.5,1),
#'     typeOfDesign = "asOF",
#'     MultipleWinners = T,
#'     Selection=F,
#'     SelectionLook = NA,
#'     SelectionMethods = NA,
#'     SelectionCriterion =NA,
#'     SelectionParmeter=NA,
#'     ImplicitSSR = 'None',
#'     nSimulation = 1000,
#'     Seed = 'Random',
#'     SummaryStat = F
#' )
#' {
#'
#'   TailType <- 'RightTail'       ##Default Right
#'   Hypothesis <- 'CommonControl' ##Default CommonControl
#'   UpdateStrategy <- F           ##Not implemented yet
#'
#'   ArmID <- 1:nArms
#'
#'   if(is.na(Arms.Name))  Arms.Name <- paste('Arm',ArmID, sep = '')
#'
#'   names(ArmID) <- names(Arms.Mean) <- names(Arms.std.dev) <- names(Arms.alloc.ratio) <- Arms.Name
#'
#'   if(test.type == 'Bonf') #Using the partly parametric function to perform Bonferroni test
#'   {
#'     Correlation <- diag(length(WI))
#'     Correlation[Correlation==0] = NA
#'   }else if(test.type == 'Sidak'||test.type == 'Simes')
#'   {
#'     Correlation <- NA
#'   }
#'
#'   #info to run Simulation
#'   gmcpSimObj <- list(
#'     #test parameters
#'     'nArms'= nArms,               'nLooks'= length(info_frac),
#'     'Max_SS' = SampleSize,        'test.type' = test.type,
#'     'IntialWeights'=WI,                      'G' = G,
#'     'Correlation' = Correlation,  'Hypothesis' = Hypothesis,
#'     'TailType' = TailType,         'alpha' = alpha,
#'
#'     #Boundary
#'     'InfoFrac' = info_frac,        'EffCutOff' = NA,
#'     'MultipleWinners'=MultipleWinners,
#'
#'     #Response Generation
#'     'ArmID' = ArmID,               'Arms.Mean' = Arms.Mean,
#'     'Arms.std.dev' = Arms.std.dev, 'Arms.alloc.ratio' = Arms.alloc.ratio,
#'     'Arms.alloc.ratio' = Arms.alloc.ratio,
#'
#'     #Selection
#'     'Selection' = Selection,               'SelectionLook' = SelectionLook,
#'     'SelectionMethods' = SelectionMethods, 'SelectionCriterion' = SelectionCriterion,
#'     'SelectionParmeter' = SelectionParmeter,
#'
#'     #Simulation Parameters
#'     'nSimulation' = nSimulation, 'Seed' = Seed,
#'     'SummaryStat' = SummaryStat,
#'
#'     #SSR
#'     'ImplicitSSR' = ImplicitSSR,
#'
#'     #Update Strategy
#'     'UpdateStrategy' =UpdateStrategy
#'
#'   )
#'
#'   mnMAMSMEP_sim(gmcpSimObj)
#'
#' }
#'
#'
