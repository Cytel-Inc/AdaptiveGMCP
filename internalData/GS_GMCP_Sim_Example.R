# GS_GMCP_Sim_Example.R
# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

############################################################################
#####   Script to Run Examples for group sequential GMCP Tests  ############
############################################################################

library(AdaptGMCP)

# Reference paper: "A Comparison of Two Methods for Adaptive Multi-Arm Two-Stage
# Design", May 2025, Mehta & Kappler
# File: Stagewise vs Cumulative MAMS.pdf

# EXAMPLE 4 #############################
# Comparing power with East
# Setting: 4 endpoints and 1 stage (fixed sample design) - parallel gatekeeping
#          Test type: Bonferroni

num_arms <- 2 # one treatment arm versus a common control
ep_type <- list(EP1 = "Continuous", EP2 = "Continuous",
                EP3 = "Continuous", EP4 = "Continuous") # 4 continuous endpoints
arm_means <- list(EP1 = c(2.3, 2.6), EP2 = c(-4.5, -2.5),
                  EP3 = c(-10, -6.5), EP4 = c(-1.2, -0.4))
arm_sd <- list(EP1 = c(1.095445115, 1.095445115),
               EP2 = c(6.480740698, 6.480740698),
               EP3 = c(12.04159458, 12.04159458),
               EP4 = c(2.828427125, 2.828427125))
alloc_ratio <- c(1, 1) # balanced design

ep_corr <- matrix(c(1, 0.507092553, 0.51550667, 0.516397779,
                    0.507092553, 1, 0.486939438, 0.507356595,
                    0.51550667, 0.486939438, 1, 0.499137187,
                    0.516397779, 0.507356595, 0.499137187, 1),
                  byrow = T, nrow = 4)

# Graph for parallel gatekeeping:
# First family: EP1 and EP2
# Second family: EP3 and EP4
# Test second family only if at least one hypothesis from the first family
# is rejected.
weights <- c(1/2, 1/2, 0, 0)
trans_mat <- matrix(c(0, 0, 1/2, 1/2,
                      0, 0, 1/2, 1/2,
                      0, 0, 0, 1,
                      0, 0, 1, 0), byrow = T, nrow = 4)

nSims <- 25000 # 10 # 100 # 1000 #

# Calling the simulation function
out4 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 250,
                   TestStatCont = "t-equal",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 4, lEpType = ep_type,
                   Arms.Mean = arm_means, Arms.std.dev = arm_sd,
                   CommonStdDev = F, Arms.alloc.ratio = alloc_ratio,
                   EP.Corr = ep_corr, WI = weights, G = trans_mat,
                   test.type = "Bonf", info_frac = 1,
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out4

# AdaptGMCP output:
# Global power = 0.68884, CI = (0.6831, 0.69458)
# Conjunctive Power = 0.25272, CI = (0.24733, 0.25811)
# Disjunctive Power = 0.68884, CI = (0.6831, 0.69458)
# East output:
# Global power = 0.68884
# Conjunctive Power = 0.20519
# Disjunctive Power = 0.68884
# Remark: Global / disjunctive powers match exactly. However, conjunctive
# powers differ.

# EXAMPLE 3 #############################
# Comparing power with East
# Setting: 1 endpoint and 1 stage (fixed sample design)
#          Test type: Dunnett Single Step

num_arms <- 4 # three treatment arms versus a common control
ep_type <- list(EP1 = "Continuous") # One continuous endpoint
arm_means <- list(EP1 = c(0, 0.45, 1.5, 3))
arm_sd <- list(EP1 = c(5, 5, 5, 5))
alloc_ratio <- c(1, 1, 1, 1) # balanced design

weights <- c(1/3, 1/3, 1/3)
trans_mat <- matrix(c(0, 1/2, 1/2,
                      1/2, 0, 1/2,
                      1/2, 1/2, 0), byrow = T, nrow = 3)

nSims <- 25000 # 10 # 100 # 1000 #

# Calling the simulation function
out3 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 200,
                   TestStatCont = "t-equal",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Mean = arm_means, Arms.std.dev = arm_sd,
                   Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Dunnett", info_frac = 1,
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out3

# AdaptGMCP output:
# Global power = 0.74844, CI = (0.74306, 0.75382)
# Conjunctive Power = 0.04404, CI = (0.0415, 0.04658)
# Disjunctive Power = 0.74844, CI = (0.74306, 0.75382)
# East output:
# Global power = 0.75138
# Conjunctive Power = 0.01948
# Disjunctive Power = 0.75138
# Remark: While global and disjunctive powers are close between AdaptGMCP and
# East, the conjunctive powers do not match.


# EXAMPLE 2 #############################
# Comparing FWER with East
# Setting: 1 endpoint and 1 stage (fixed sample design)
#          Test type: Weighted Bonferroni

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Continuous") # One continuous endpoint
arm_means <- list(EP1 = c(5, 5, 5, 5, 5))
arm_sd <- list(EP1 = c(10, 10, 10, 10, 10))
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- c(0.4, 0.4, 0.1, 0.1)
trans_mat <- matrix(c(0, 1/3, 1/3, 1/3,
                      1/3, 0, 1/3, 1/3,
                      1/3, 1/3, 0, 1/3,
                      1/3, 1/3, 1/3, 0), byrow = T, nrow = 4)

nSims <- 25000 # 10 # 100 # 1000 #

# Calling the simulation function
out2 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                   TestStatCont = "t-equal",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Mean = arm_means, Arms.std.dev = arm_sd,
                   Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Bonf", info_frac = 1,
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out2

# FWER = 0.0212, CI = (0.018, 0.024)


# EXAMPLE 1 #############################
# Simulations >>>>>>>>>>>>>>>>>>>>>>>>>>
# From page 11 of the paper
# 5-arm trial (including 1 control), balanced randomization
# Data generated from Normal distribution with sigma=1
# Total sample size for the trial = 1000 (200 per arm)
# 2 look design, look 1 after 500 subjects (100 per arm) observed
# alpha=0.025, LDOF spending function
# Trial terminates at look 1 if any Hj is rejected by closed testing.
# Otherwise arms are dropped / selected according to one of pre-specified 4 rules.
# Under three of these rules (conservative, normal, and agressive), there is a
# varying selection threshold. If all Wald stats at look 1 fall below this
# threshold, all treatment arms are dropped and the trial stops for futility.
# The 4th rule selects the best dose at the interim look. So, three of the four
# treatment arms are always dropped and the trial always continues to the second
# look under that rule. The trial cannot stop early for futility in that case.
# At stage 2, the remaining 500 subjects are assigned in equal proportion to
# the selected arms. This constitutes implicit SS increase if arms were dropped
# at stage 1.

# Setting the input
# First testing with the p-value combination (stagewise MAMS) method
meth <- "CER" # CombPValue" # Param Method
alp <- 0.025
SS <- 1000
stat <- "t-equal" # Param TestStatCont
# Ignore TestStatBin

fwer_ctrl <- "None" # CombinationTest"
# The documentation says that FWERControl is applicable for CER method only. But
# it takes one of the 2 values: "CombinationTest" which uses incremental test
# stats and and "None", which uses cumulative test stats. This sounds redundant
# and potentially problematic since the parameter Method already specifies whether
# you are using the stagewise or cumulative MAMS. So, why is CombinationTest
# required again?
# TODO: Investigate this.

arms <- 5 # including control
endpts <- 1 # only 1 endpoint in this case
ep_type <- list(EP1="Continuous")

# SCENARIO 1: is FWER preserved under global null hypo?
arm_mean <- list(EP1=c(0, 0, 0, 0, 0))

# # SCENARIO 2: one efficacious treatment
# arm_mean <- list(EP1=c(0, 0, 0, 0, 0.25))

arm_sd <- list(EP1=c(1, 1, 1, 1, 1)) # Param Arms.std.dev

# TODO: Investigate how Arms.std.dev and CommonStdDev work together. How do we
# specify a single std dev value of CommonStdDev is TRUE?
# It appears that in that case, you are supposed to create a list of only 1
# element with the common sd value for Arms.std.dev.
comm_sd <- T

# Ignore Arms.Prop

arm_alloc <- c(1, 1, 1, 1, 1)

# # TODO: Investigate if EP.Corr param specifies correlation between endpoints
# # when multiple endpoints are specified or it specifies correlation between
# # the different test stats like in case of analysis with stagewise MAMS.
# # Assuming that it is for multiple endpoints as the help page says, ignore it
# # for now since we have only one endpoint in this example.
# # corr <- matrix(c(1, 0.5, 0.5, 1), byrow = T, nrow = 2) # correlation matrix
# corr <- matrix(c(1, 0.5, 0.5, 0.5,
#                  0.5, 1, 0.5, 0.5,
#                  0.5, 0.5, 1, 0.5,
#                  0.5, 0.5, 0.5, 1), byrow = T, nrow = 4)

wi <- rep(1/4,4) # Weights for the 4 hypo in the initial graph
# Transition matrix for the initial graph
g <- matrix(c(0, 1/3, 1/3, 1/3,
              1/3, 0, 1/3, 1/3,
              1/3, 1/3, 0, 1/3,
              1/3, 1/3, 1/3, 0), byrow = T, nrow = 4)

# The paper does not specify what MCP test is to be used. Trying with Bonferroni
tt <- "Parametric" # Bonf" # test.type

t <- c(0.5, 1) # info fractions for looks 1 and 2
des <- "asOF"
# Parameters deltaWT, deltaPT1, gammaA, userAlphaSpending ignored

# Param MultipleWinners
mult_win <- F # Stop the trial if any hypo is rejected at look 1

sel <- T # Selection required at look 1
sel_lk <- 1 # Look at which selection happens
sel_ep <- 1 # Which endpoint to select (there is only 1 in this example)
sel_scale <- "teststat" # Param SelectionScale: arms with test stat exceeding
                        # the given threshold are selected

# RULE 1: Conservative drop-the-loser rule: Drop treatment j if zj < −0.6745
# At look 2, remaining 500 subjects are allocated in equal proportion to
# selected arms.
sel_crit <- "threshold" # Param SelectionCriterion
sel_param <- -0.6745 # Param SelectionParmeter
# TODO: Investigate how the package understand whether teststat < threshold OR
# teststat > threshold is to be used as there is no obvious input parameter
# indicating the comparison type, which should make a difference for left-tailed
# vs right-tailed tests.

keep <- T # whether to keep all the associated hypothesis for the selected arms
# TODO: Investigate what happens if this is set to FALSE

# Param ImplicitSSR
impl_SSR <- "Selection" # re-allocate samples only from de-selected arms to
# available arms
# TODO: It is not clear how tthe ImplicitSSR is implemented since the help
# message is confusing. Investigate.

nSim <- 1000
# Ignoring nSimulation_Stage2 as it is not applicable.
# Ignore Seed so that it is left as random.

summ_stat <- F # Summary Stat file not required

graph <- T

East_stat <- NULL # East summary stat file not required
parallel <- F

# Ignore UpdateStrategy

# Calling function...
# Without EP.Corr param
simMAMSMEP(Method = meth, alpha = alp, SampleSize = SS, TestStatCont = stat,
           FWERControl = fwer_ctrl, nArms = arms, nEps = endpts, lEpType = ep_type,
           Arms.Mean = arm_mean, Arms.std.dev = arm_sd, CommonStdDev = comm_sd,
           Arms.alloc.ratio = arm_alloc, # EP.Corr = corr,
           WI = wi, G=g, test.type = tt,
           info_frac = t, typeOfDesign = des, MultipleWinners = mult_win,
           Selection = sel, SelectionLook = sel_lk, SelectEndPoint = sel_ep,
           SelectionScale = sel_scale, SelectionCriterion = sel_crit,
           SelectionParmeter = sel_param, KeepAssosiatedEps = keep,
           ImplicitSSR = impl_SSR, nSimulation = nSim, SummaryStat = summ_stat,
           plotGraphs = graph, EastSumStat = East_stat, Parallel = parallel)

# With EP.Corr param
# simMAMSMEP(Method = meth, alpha = alp, SampleSize = SS, TestStatCont = stat,
#            FWERControl = fwer_ctrl, nArms = arms, nEps = endpts, lEpType = ep_type,
#            Arms.Mean = arm_mean, Arms.std.dev = arm_sd, CommonStdDev = comm_sd,
#            Arms.alloc.ratio = arm_alloc, EP.Corr = corr,
#            WI = wi, G=g, test.type = tt,
#            info_frac = t, typeOfDesign = des, MultipleWinners = mult_win,
#            Selection = sel, SelectionLook = sel_lk, SelectEndPoint = sel_ep,
#            SelectionScale = sel_scale, SelectionCriterion = sel_crit,
#            SelectionParmeter = sel_param, KeepAssosiatedEps = keep,
#            ImplicitSSR = impl_SSR, nSimulation = nSim, SummaryStat = summ_stat,
#            plotGraphs = graph, EastSumStat = East_stat, Parallel = parallel)

#########################################################

# TEMP CODE ##########################


