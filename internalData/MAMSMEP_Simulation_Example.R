
library(AdaptGMCP)

timeNow <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
sink(paste0("Console_output_", timeNow, ".txt"), append = TRUE, type = 'output')
sink(paste0("Console_output_", timeNow, ".txt"), append = TRUE, type = 'message')

cat(readChar(rstudioapi::getSourceEditorContext()$path, # Writing currently opened R script to file
             file.info(rstudioapi::getSourceEditorContext()$path)$size))

########################## Inputs ##############################
#------------------Test Method----------------------------------
#options- 'CombPValue': Combining p-values, 'CER': Conditional Error
Method <- 'CER' #CombPValue/CER

#-----------------------Total Sample Size------------------------
SampleSize <- 324
#-----------------------Design Alpha-----------------------------
alpha <- 0.025

#----------------------Test Statistics---------------------------
#options : 'z' :- Z Statistics, 't-equal' :- t Statistics with equal variance, 't-unequal' :- t Statistics with unequal variance
TestStat <- 't-equal'

#---------FWER Control Method(Applicable for CER Method only)------
#options: 'CombinationTest' : combined incremental p-values for stage-2 testing, 'None': Cumulative classical test statistics.
FWERControl <- 'None'

#-----------------------Number of Arms--------------------------
nArms <- 3

#-------------------Number of Endpoints-------------------------
nEps  <- 2
#-------------Arm-Wise mean for each endpoints-----------------
Arms.Mean <- list('EP1' = c(0,0.4,0.3),
                 'EP2' = c(0.1,0.45,0.25))

#-------------Arm-Wise planned Std. Dev. for each endpoints-----
Arms.std.dev <- list('EP1' = c(1.1,1.2,1.3),
                     'EP2' = c(1.5,1.6,1.7))

#-------------Arm-Wise planned allocation Ratio----------------
Arms.alloc.ratio <- c(1,1,1)

#-------------Correlation between End points -------------------
EP.Corr <- matrix(c(1,0.5,
                   0.5,1),
                 nrow = 2)
#----------------Initial Weights---------------------
WI <-  c(0.5,0.5,0,0)

#---------------Transition Matrix--------------------
G <- matrix(c(0,0.5,0.5,0,
             0.5,0,0,0.5,
             0,1,0,0,
             1,0,0,0),
           nrow = nEps*(nArms-1), byrow = T)

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
SelectionScale <- 'pvalue'

#-----------------Selection Criteria----------------------
#Options- 'best': best r, 'threshold': threshold for selection, 'epsilon': for epsilon neighborhood
SelectionCriterion <- 'best'

#-----------------Selection Criteria Parameter----------------------
# r for best, threshold value for threshold and epsilon distance for criteria 'epsilon'
SelectionParmeter <- 1

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

#-----------------------Run Simulation--------------------------------------
#Please uncomment the following code to run the simulation(short-cut to uncomment 1.Select the lines, 2.ctr+shift+c)
#
out <- simMAMSMEP(alpha = alpha, SampleSize = SampleSize,nArms = nArms,nEps = nEps, TestStat=TestStat, FWERControl = FWERControl,
           Arms.Mean = Arms.Mean, Arms.std.dev = Arms.std.dev, Arms.alloc.ratio = Arms.alloc.ratio,
           EP.Corr = EP.Corr,WI = WI, G = G, test.type = test.type,info_frac = info_frac,
           typeOfDesign = typeOfDesign, MultipleWinners = MultipleWinners,
           Selection = Selection,SelectionLook = SelectionLook,SelectEndPoint = SelectEndPoint,SelectionScale = SelectionScale,
           SelectionCriterion = SelectionCriterion, SelectionParmeter = SelectionParmeter, KeepAssosiatedEps = KeepAssosiatedEps,
           ImplicitSSR = ImplicitSSR, nSimulation = nSimulation, Seed = Seed, SummaryStat = SummaryStat,
           Method = Method, plotGraphs = plotGraphs)
out

# sink()
closeAllConnections() # Close connection to log file
save(out, file = paste0("Output_simMAMSMEP_",timeNow, ".RData"))

# # Use the following code to reload the file
# load("Result_simMAMSMEP_2024-01-31_12-29-48.RData")
