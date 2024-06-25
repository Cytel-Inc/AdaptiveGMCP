
#######################################################################
###############  Script to Run Adapt GMCP Simulations #################
#######################################################################

library(AdaptGMCP)
########################## Inputs ##############################
#------------------Test Method----------------------------------
#options- 'CombPValue': Combining p-values, 'CER': Conditional Error
Method <- 'CER'

#-----------------------Total Sample Size------------------------
SampleSize <- 400

#-----------------------Design Alpha-----------------------------
alpha <- 0.025

#-----------------------Number of Arms--------------------------
#Note: Including Control Arm
nArms <- 4

#-------------------Number of Endpoints-------------------------
nEps  <- 2

#-----------------------Endpoint Type----------------------------
# List containing the endpoint types
# Options : "Continuous", "Binary"
EpType <- list("EP1" = "Continuous",
               "EP2" = "Binary")

#----------------------Test Statistics(Continuous)---------------------------
#options : "z" :- Z Statistics, "t-equal" :- t Statistics with equal variance, "t-unequal" :- t Statistics with unequal variance
TestStatCon <- "t-equal"

#----------------------Test Statistics(Binary)---------------------------
#options : "Pooled" :- Z Statistics with equal variances, "UnPooled" :- Z Statistics with un-equal variances
TestStatBin <- "UnPooled"

#---------FWER Control Method(Applicable for CER Method only)------
#options: 'CombinationTest' : combined incremental p-values for stage-2 testing, 'None': Cumulative classical test statistics.
FWERControl <- "None"


#-------------Arm-Wise mean for each endpoints-----------------
# List containing the arm-wise(continuous) means
#Note: The first input is for control arm and the rest are for the treatments
Arms.Mean <- list('EP1' = c(0, 0.4, 0.4, 0.4),
                  'EP2' = NA)

#-------------Arm-Wise planned Std. Dev. for each endpoints-----
# List containing the arm-wise(continuous) standard deviations
# Note: The first input is for control arm and the rest are for the treatments
Arms.std.dev <- list('EP1' = c(1.1, 1.2, 1.3, 1.4),
                     'EP2' = NA)

#--------Use Common standard deviations for boundary computation-----
#TRUE = the treatment standard deviations assumed to be same as the control for boundary computations for continuous endpoints
#FALSE = the treatment standard deviations assumed to be same as given in Arms.std.dev.
CommonStdDev <- FALSE


#-------------Arm-Wise proportion for each endpoints-----------------
# List containing the arm-wise(Binary) proportions
# Note: The first input is for control arm and the rest are for the treatments
Arms.Prop <- list('EP1'=NA,
                  'EP2' =c(0.1, 0.4, 0.4, 0.4))

#-------------Arm-Wise planned allocation Ratio----------------
# Vector containing the arm-wise allocation ratios
# Note: The first input is for control arm and the rest are for the treatments
Arms.alloc.ratio <- c(1, 1, 1, 1)

#-------------Correlation between End points -------------------
# Note: For composite endpoints this correlation is assumed to be normal scale correlations
EP.Corr <- matrix(c(1, 0.5,
                    0.5, 1),
                  nrow = nEps)

#----------------Initial Weights---------------------
#Note: Hypotheses will follow the order of Endpoints and Treatments as given in 'Arms.Mean' and 'Arms.std.dev' inputs
#e.g.: If 'Arms.Mean' are given in the format list('EP1'=c(ctr_mean, trt1_mean, trt2_mean), 'EP2'=c(ctr_mean, trt1_mean, trt2_mean)) then the hypothesis will be
# H1 = (Trt1 vs Ctr for EP1), H2 = (Trt2 vs Ctr for EP1), H3 = (Trt1 vs Ctr for EP2), H4 = (Trt2 vs Ctr for EP2)
#The initial weights and the transition matrix will follow the order accordingly as (H1,H2,H3,H4)

WI <-  c(1/3, 1/3, 1/3, 0, 0, 0)

#---------------Transition Matrix--------------------
G <- matrix(c(0,0,0,     1,0,0,
              0,0,0,     0,1,0,
              0,0,0,     0,0,1,
              0,1/2,1/2, 0,0,0,
              1/2,0,1/2, 0,0,0,
              1/2,1/2,0, 0,0,0),
            nrow = nEps*(nArms-1), byrow = TRUE)



#-------------------Test Type------------------------
# Testing Procedure
#Available tests for Combining P-values Method :- 'Bonf','Sidak','Simes','Dunnett','Partly-Parametric'
#Available tests for CER Method :- 'Parametric','Non-Parametric','Partly-Parametric'
test.type <- 'Partly-Parametric'

#----------------Information Fraction----------------
#The number of lookes is same as the length of info_frac(for FSD info_frac = 1)
info_frac <-  c(1/2,1)

#-------------------Boundary Type-------------------
# O'Brien & Fleming = 'OF', Pocock = 'P', Wang & Tsiatis Delta class = 'WT'
# Pampallona & Tsiatis = 'PT', Haybittle & Peto = 'HP',
# Optimum design within Wang & Tsiatis class = 'WToptimum'
# O'Brien & Fleming type alpha spending = 'asOF', Pocock type alpha spending = 'asP'
# Kim & DeMets alpha spending = 'asKD', Hwang, Shi & DeCani alpha spending = 'asHSD'
# no early efficacy stop = 'noEarlyEfficacy'
# default is "asOF"

typeOfDesign <- "asOF"

#-----------------Multiple Winner Option--------------
MultipleWinners <- TRUE

#-----------------Treatment Selection Choice--------------
Selection <- TRUE

#-----------------Selection based on look-------------
SelectionLook <- 1

#-----------------Selection Endpoint----------------------
#options- '1': Endpoint 1, '2':Endpoint 2, 'overall': overall
SelectEndPoint <- 1

#-----------------Selection Based on----------------------
#Options- 'delta': delta, 'teststat': Test Statistics, 'stderror' : Standard Error of the test stat,  'pvalue': p-value(un-adj)
SelectionScale <- 'teststat'

#-----------------Selection Criteria----------------------
#Options- 'best': best r, 'threshold': threshold for selection, 'epsilon': for epsilon neighborhood
SelectionCriterion <- 'threshold'

#-----------------Selection Criteria Parameter----------------------
# r for best, threshold value for threshold and epsilon distance for criteria 'epsilon'
SelectionParmeter <- 0.6745

#----------------Keep associated Hypothesis-------------------------
#TRUE: keep all the associated hypothesis for the selected arms, FALSE otherwise
KeepAssosiatedEps <- TRUE

#-------------------Implicit SSR-------------------------------------
#'Selection': re-allocate samples only from de-selected arms to available arms,
#'All': Allocate all the planned samples(for the look) to the available arms,
#'None': No Re-allocation
ImplicitSSR <- 'Selection'

#-------------------Number of Simulations-------------------------------------
nSimulation <- 1000

#-------------------Simulation Seed-------------------------------------
Seed <- 100

#------------------Print Summary Statistics file------------------------
SummaryStat <- FALSE

#-------------------Plot Initial Graph------------------------
plotGraphs <- TRUE

#--------------Parallel Logical; TRUE: Parallel computations--------
Parallel <- TRUE

#-----------------------Run Simulation--------------------------------------
# Uncomment the following code to run simulation;
# short cut: 1) select all the following lines 2) ctrl+shift+c]
#
# out <- simMAMSMEP(
#   alpha = alpha, SampleSize = SampleSize, nArms = nArms, nEps = nEps,lEpType=EpType,
#   TestStatCon = TestStatCon, TestStatBin = TestStatBin, FWERControl = FWERControl,
#   Arms.Mean = Arms.Mean, Arms.std.dev = Arms.std.dev, CommonStdDev = CommonStdDev, Arms.Prop = Arms.Prop, Arms.alloc.ratio = Arms.alloc.ratio,
#   EP.Corr = EP.Corr, WI = WI, G = G, test.type = test.type, info_frac = info_frac,
#   typeOfDesign = typeOfDesign, MultipleWinners = MultipleWinners,
#   Selection = Selection, SelectionLook = SelectionLook, SelectEndPoint = SelectEndPoint, SelectionScale = SelectionScale,
#   SelectionCriterion = SelectionCriterion, SelectionParmeter = SelectionParmeter, KeepAssosiatedEps = KeepAssosiatedEps,
#   ImplicitSSR = ImplicitSSR, nSimulation = nSimulation, Seed = Seed, SummaryStat = SummaryStat,
#   Method = Method, plotGraphs = plotGraphs, Parallel = Parallel
# )
# out


