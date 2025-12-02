
#######################################################################
###############   Script to Run CER Analysis examples #################
#######################################################################

library(AdaptGMCP)

########################## Inputs ##############################
#-----------------------Total Sample Size------------------------
SampleSize <- 210

#-----------------------Design Alpha------------------------
alpha <- 0.025

#-----------------------Number of Arms------------------------
#Note: Including Control Arm
nArms <- 3

#-------------------Number of Endpoints------------------------
nEps  <- 1

#------------------list with endpoint types--------------------
# List containing the endpoint types
# Options : "Continuous", "Binary"
EpType <- list("EP1" = "Continuous")


#-------------Arm-Wise planned Std. Dev. for each endpoints-----
# List containing the arm-wise(continuous) standard deviations
# The First input is for the control arm
# Not required for test.type = 'Non-Parametric' or endpoint type as "Binary"
sigma <- list(
  "EP1" = c(2, 2, 2) # c(2, 5, 8) # c(1, 5, 8) # c(1, 2, 3) # c(1, 1, 1) # c(0.01, 0.5, 0.1) #
  )

#--------Use Common standard deviations for boundary computation-----
#TRUE = the treatment standard deviations assumed to be same as the control for boundary computations for continuous endpoints
#FALSE = the treatment standard deviations assumed to be same as given in Arms.std.dev.
CommonStdDev <- FALSE

#-------------proportion for control arm-------------------
# List containing the Control Arm(Binary) Proportions
# Not required for test.type = 'Non-Parametric' or endpoint type as "Continuous"
prop.ctr <- list("EP1" = NA,
                 "EP2" = 0.2)

#-----------------------Arm-Wise Allocation Ratio------------------------
# Vector of length = nArms
# Not required for test.type = 'Non-Parametric'
allocRatio <- c(1,1,1)

#----------------Initial Weights---------------------
#The number of hypothesis=(nEps*(nArms-1)) is same as the length of WI
WI <-  c(1/2, 1/2)

#---------------Transition Matrix--------------------
#The dimension of the matrix is (Number of Hypothesis x Number of Hypothesis)
G <- matrix(c(0,1,
              1,0),
      nrow = nEps*(nArms-1), byrow = TRUE)

#----------------Information Fraction----------------
#The number of lookes is same as the length of info_frac(for FSD info_frac = 1)
info_frac <- c(0.5,1)

#-------------------Boundary Type-------------------
# O'Brien & Fleming = 'OF', Pocock = 'P', Wang & Tsiatis Delta class = 'WT'
# Pampallona & Tsiatis = 'PT', Haybittle & Peto = 'HP',
# Optimum design within Wang & Tsiatis class = 'WToptimum'
# O'Brien & Fleming type alpha spending = 'asOF', Pocock type alpha spending = 'asP'
# Kim & DeMets alpha spending = 'asKD', Hwang, Shi & DeCani alpha spending = 'asHSD'
# no early efficacy stop = 'noEarlyEfficacy'
# default is "asOF"

typeOfDesign <- 'asOF'

#-------------------Test Type------------------------
# Testing Procedure
# options = 'Parametric', 'Partly-Parametric', 'Non-Parametric'
# Note: For parametric tests the the number of endpoints must be 1, whereas for partly parametric case it has to be greater than 1
test.type <- 'Parametric'


#-------------------Adaptation Flag------------------------
# AdaptStage2 : TRUE = Adaptation options will be given after stage-1
AdaptStage2 <- TRUE

#-------------------Plot Intermediate Graphs------------------------
plotGraphs <- TRUE

#--------------Run Analysis--------------------------
###For Interim-look inputs follow the R console#####
#Please uncomment the following code to run the Analysis(short-cut to uncomment 1.Select the lines, 2.ctr+shift+c)
#
#

#sink("Experimental/output.txt")

adaptGMCP_CER(nArms = nArms, nEps = nEps, EpType = EpType,
              sigma = sigma, CommonStdDev = CommonStdDev,
              prop.ctr = prop.ctr, allocRatio = allocRatio,
              SampleSize = SampleSize, alpha = alpha, WI = WI,
              G = G, info_frac = info_frac, typeOfDesign = typeOfDesign,
              test.type = test.type, AdaptStage2 = AdaptStage2,
              plotGraphs = plotGraphs)
