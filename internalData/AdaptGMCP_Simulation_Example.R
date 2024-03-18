#######################################################################
#########   Script to Run adaptive gmcp Simulation examples ###########
#######################################################################

#-------------------Sample Size---------------------
SampleSize <- 500

#-----------------------Alpha------------------------
alpha <- 0.025

#----------------Information Fraction----------------
#The number of lookes is same as the length of info_frac
#For one look specify as 1
info_frac <- c(1/2,1)

#-------------------Boundary Type-------------------
# O'Brien & Fleming = 'OF', Pocock = 'P', Wang & Tsiatis Delta class = 'WT'
# Pampallona & Tsiatis = 'PT', Haybittle & Peto = 'HP',
# Optimum design within Wang & Tsiatis class = 'WToptimum'
# O'Brien & Fleming type alpha spending = 'asOF', Pocock type alpha spending = 'asP'
# Kim & DeMets alpha spending = 'asKD', Hwang, Shi & DeCani alpha spending = 'asHSD'
# no early efficacy stop = 'noEarlyEfficacy'
# default is "asOF"

typeOfDesign <- 'asOF'

#------------------Response Generation Parameters------------------
#Number of Arms(Including control)
nArms <- 5

#Arm-Wise Mean
Arms.Mean <- c(0,0.3,0.3,0.1,0.0)

#Arm-Wise Standard Deviations
Arms.std.dev <- c(1,1,1,1,1)

#Arm-Wise Allocation Ratio
Arms.alloc.ratio <- c(1,1,1,1,1)

#----------------Initial Weights---------------------
#The number of hypothesis(k) is same as the length of WI
WI <- c(1/2,1/2,0,0)

#---------------Transition Matrix--------------------
#The dimension of the matrix is (Number of Hypothesis x Number of Hypothesis)
G <- matrix(c(0,0.75,0.25,0,
              0.75,0,0,0.25,
              0,1,0,0,
              1,0,0,0),
            nrow = 4, byrow = TRUE)

#------------Plot the Testing Strategy---------------
#Group of hypothesis
hGroup <- c('P','P','S','S')

#Coordinates of the nodes(Hypothesis)
cordinates <- list(c(-1,1),c(1,1),c(-1,-1),c(1,-1))

#Plot
gmcpPlot(WI = WI, G = G, hGroup = hGroup, cordinates = cordinates)

#-------------------Test Type------------------------
# Testing Procedure
# Bonferroni = 'Bonf' , Sidak = 'Sidak', Simes ='Simes'
# Dunnett = 'Dunnett', Partly Parametric = 'Partly-Parametric'

test.type <- 'Partly-Parametric'


#---------------Correlation--------------------------
#Need to modify only for Dunnett and Partly Parametric test else will be ignored inside the code.
#The dimension of the matrix is (Number of Hypothesis x Number of Hypothesis)
Correlation <- matrix(c(1,0.5,NA,NA,
                        0.5,1,NA,NA,
                        NA,NA,1,0.5,
                        NA,NA,0.5,1),
                      byrow=TRUE, nrow= 4)


#--------------Selection-----------------------------
#Selection Choice; True if selection required at interim.
Selection <- TRUE

#SelectionLook: ScalerVector to specify the selection look/looks
SelectionLook <- 1

#SelectionMethods: Scale parameter on which selection will be based on. 'pvalue': p-value(un-adj) based selection (Currently available)
SelectionMethods <- 'pvalue'

#SelectionCriterion: best' = best r, 'threshold' = threshold for selection
SelectionCriterion = 'best'

#SelectionParmeter: r for best, threshold value for threshold
SelectionParmeter <- 1

#--------------------Multiple Winner & Implicit SSR---------------------------

#MultipleWinners = TRUE if multiple winners FALSE otherwise.
MultipleWinners <- TRUE

#ImplicitSSR:
#'Selection' if re-allocate samples only from de-selected arms to available arms
#'All': Allocate all the planned samples(for the look) to the available arms
#'None': No Re-allocation

ImplicitSSR <- 'All'

#---------------------------Simulation Parameters-------------------------------
#Number of simulation
nSimulation <- 1000

#Seed
Seed <- 100

#Option to see simulation-wise data
SummaryStat <- FALSE

#-----------------------Run Simulation--------------------------------------
#Please uncomment the following code to run the simulation(short-cut to uncomment 1.Select the lines, 2.ctr+shift+c)

# Results <- adaptGMCP_SIM(
#   SampleSize=SampleSize,alpha=alpha, nArms=nArms,Arms.Mean=Arms.Mean,Arms.std.dev=Arms.std.dev
#   ,Arms.alloc.ratio=Arms.alloc.ratio,WI=WI,G=G,test.type=test.type, Correlation = Correlation,
#   info_frac=info_frac,typeOfDesign = typeOfDesign, MultipleWinners=MultipleWinners,
#   Selection=Selection,SelectionLook=SelectionLook, SelectionMethods=SelectionMethods, SelectionCriterion=SelectionCriterion,SelectionParmeter=SelectionParmeter
#   ,ImplicitSSR=ImplicitSSR,
#   nSimulation=nSimulation,Seed=Seed,SummaryStat=SummaryStat
# )
# Results
