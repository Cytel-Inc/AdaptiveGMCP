# File: GS_GMCP_Sim_Bin_Example.R
# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

# Script for testing simMAMSMEP() for fixed sample designs witth a binomial endpoint

library(AdaptGMCP)

#########################################################################
############################ TESTS FOR POWER ############################
############################ SINGLE ENDPOINT ############################
#########################################################################

# TEST 5-POW #########################
# Setting: 1 endpoint and 2 stages
#          Test type: Weighted Simes
#          Early stopping with OF efficacy boundary
#          Multiple winners not allowed (East only has this option)

num_arms <- 4 # three treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = c(0.1, rep(0.22, 3)))
alloc_ratio <- c(1, 1, 1, 1) # balanced design

weights <- c(0.5, 0.3, 0.2)
trans_mat <- matrix(c(0, 2/3, 1/3,
                      0, 0, 1,
                      1/2, 1/2, 0), byrow = T, nrow = 3)

t <- c(4/7, 1)

nSims <- 25000 # 10 # 100 # 1000 #

# Calling the simulation function
out5.pow <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 700,
                   TestStatBin = "Pooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Simes", info_frac = t,
                   typeOfDesign = "asOF",
                   MultipleWinners = F, Selection = T, SelectionScale = "delta",
                   SelectionCriterion = "best", SelectionParmeter = 2,
                   ImplicitSSR = "All",
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out5.pow

# Results:
# AdaptGMCP output: Global power = Disjunctive power = 0.95588,
# CI = (0.95333, 0.95843)
# Conjunctive power = 0.12308, CI = (0.11901, 0.12715)
# East output:

# TEST 4-POW #########################
# Setting: 1 endpoint and 2 stages
#          Test type: Simes
#          Early stopping with OF efficacy boundary
#          Multiple winners not allowed (East only has this option)

num_arms <- 4 # three treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = c(0.1, rep(0.22, 3)))
alloc_ratio <- c(1, 1, 1, 1) # balanced design

weights <- rep(1/3, 3)
trans_mat <- matrix(c(0, 1/2, 1/2,
                      1/2, 0, 1/2,
                      1/2, 1/2, 0), byrow = T, nrow = 3)

t <- c(4/7, 1)

nSims <- 25000 # 100 # 1000 # 10 #

# Calling the simulation function
out4.pow <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 700,
                   TestStatBin = "Pooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Simes", info_frac = t,
                   typeOfDesign = "asOF",
                   MultipleWinners = F, Selection = T, SelectionScale = "delta",
                   SelectionCriterion = "best", SelectionParmeter = 2,
                   ImplicitSSR = "All",
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out4.pow

# Results:
# AdaptGMCP output:
# East output:

# TEST 3-POW #########################
# Setting: 1 endpoint and 2 stages
#          Test type: Bonferroni
#          No early stopping
#          Multiple winners not allowed (East only has this option)

num_arms <- 4 # three treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = c(0.1, rep(0.22, 3)))
alloc_ratio <- c(1, 1, 1, 1) # balanced design

weights <- rep(1/3, 3)
trans_mat <- matrix(c(0, 1/2, 1/2,
                      1/2, 0, 1/2,
                      1/2, 1/2, 0), byrow = T, nrow = 3)

t <- c(2/3, 1)

nSims <- 25000 # 100 # 1000 # 10 #

# Calling the simulation function
out3.pow <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 600,
                   TestStatBin = "UnPooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Bonf", info_frac = t,
                   typeOfDesign = "noEarlyEfficacy",
                   MultipleWinners = F, Selection = T, SelectionScale = "delta",
                   SelectionCriterion = "best", SelectionParmeter = 1,
                   ImplicitSSR = "All",
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out3.pow

# Results:
# AdaptGMCP output:
# East output:


# Example provided by Cyrus:
# For the binomial example with four treatment arms versus a common control,
# use the following design inputs. Two stage design. alpha=0.025, one sided.
# N=100/arm with interim analysis at 50/arm and LDOF spending function (i.e.
# spend 0.0015253 of the available 0.025 at stage 1).
# pi_0=0.1. pi_1=pi_2=pi_3=pi_4=0.22.
# In East MAMS this has 0.815 power. You can start with this and then add the
# weights and then add the secondary endpoints. For the secondary endpoints,
# use the same values of pi. Assume a correlation of 0.5 between the primary
# and secondary endpoints.

# TEST 2-POW #########################
# Setting: 1 endpoint and 2 stages
#          Test type: Sidak
#          Early stopping with OF efficacy boundary
#          Multiple winners not allowed (East only has this option)

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = c(0.1, rep(0.22, 4)))
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- rep(1/4, 4)
trans_mat <- matrix(c(0, 1/3, 1/3, 1/3,
                      1/3, 0, 1/3, 1/3,
                      1/3, 1/3, 0, 1/3,
                      1/3, 1/3, 1/3, 0), byrow = T, nrow = 4)

t <- c(0.5, 1)

nSims <- 25000 # 100 # 1000 # 10 #

# Calling the simulation function
out2.pow <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                   TestStatBin = "Pooled", # UnPooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Sidak", info_frac = t,
                   typeOfDesign = "asOF",
                   MultipleWinners = F, Selection = T, SelectionScale = "delta",
                   SelectionCriterion = "best", SelectionParmeter = 1,
                   ImplicitSSR = "All",
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out2.pow

# Results:
# AdaptGMCP output:
# East output:

# TEST 1-POW #########################
# Comparing power with East
# Setting: 1 endpoint and 2 stages
#          Test type: Sidak
#          No early stopping

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = c(0.1, rep(0.22, 4)))
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- rep(1/4, 4)
trans_mat <- matrix(c(0, 1/3, 1/3, 1/3,
                      1/3, 0, 1/3, 1/3,
                      1/3, 1/3, 0, 1/3,
                      1/3, 1/3, 1/3, 0), byrow = T, nrow = 4)

t <- c(0.5, 1)

nSims <- 25000 # 100 # 1000 # 10 #

# Calling the simulation function
out1.pow <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                   TestStatBin = "Pooled", # UnPooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Sidak", info_frac = t,
                   typeOfDesign = "noEarlyEfficacy",
                   MultipleWinners = F, Selection = T,
                   SelectionCriterion = "best", SelectionParmeter = 1,
                   ImplicitSSR = "All",
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out1.pow

# Results:
# AdaptGMCP output:
# East output:


########################################################################
############################ TESTS FOR FWER ############################
############################ SINGLE ENDPOINT ###########################
########################################################################

# TEST 5 #########################
# Setting: 1 endpoint and 2 stages
#          Test type: Weighted Simes
#          Early stopping with OF efficacy boundary
#          Multiple winners not allowed (East only has this option)

num_arms <- 4 # three treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = rep(0.22, 4)) # global null hypothesis
alloc_ratio <- c(1, 1, 1, 1) # balanced design

weights <- c(0.5, 0.3, 0.2)
trans_mat <- matrix(c(0, 2/3, 1/3,
                      0, 0, 1,
                      1/2, 1/2, 0), byrow = T, nrow = 3)

t <- c(4/7, 1)

nSims <- 25000 # 10 # 100 # 1000 #

# Calling the simulation function
out5 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 700,
                   TestStatBin = "Pooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Simes", info_frac = t,
                   typeOfDesign = "asOF",
                   MultipleWinners = F, Selection = T, SelectionScale = "delta",
                   SelectionCriterion = "best", SelectionParmeter = 2,
                   ImplicitSSR = "All",
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out5

# Results:
# AdaptGMCP output:
# East output:

# TEST 4 #########################
# Setting: 1 endpoint and 2 stages
#          Test type: Simes
#          Early stopping with OF efficacy boundary
#          Multiple winners not allowed (East only has this option)

num_arms <- 4 # three treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = rep(0.22, 4)) # global null hypothesis
alloc_ratio <- c(1, 1, 1, 1) # balanced design

weights <- rep(1/3, 3)
trans_mat <- matrix(c(0, 1/2, 1/2,
                      1/2, 0, 1/2,
                      1/2, 1/2, 0), byrow = T, nrow = 3)

t <- c(4/7, 1)

nSims <- 25000 # 100 # 1000 # 10 #

# Calling the simulation function
out4 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 700,
                   TestStatBin = "Pooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Simes", info_frac = t,
                   typeOfDesign = "asOF",
                   MultipleWinners = F, Selection = T, SelectionScale = "delta",
                   SelectionCriterion = "best", SelectionParmeter = 2,
                   ImplicitSSR = "All",
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out4

# Results:
# AdaptGMCP output:
# FWER = 0.019, CI = (0.01731, 0.02069)
# East output:
# FWER = 0.01974

# TEST 3 #########################
# Setting: 1 endpoint and 2 stages
#          Test type: Bonferroni
#          Early stopping with OF efficacy boundary
#          Multiple winners not allowed (East only has this option)

num_arms <- 4 # three treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = rep(0.22, 4)) # global null hypothesis
alloc_ratio <- c(1, 1, 1, 1) # balanced design

weights <- rep(1/3, 3)
trans_mat <- matrix(c(0, 1/2, 1/2,
                      1/2, 0, 1/2,
                      1/2, 1/2, 0), byrow = T, nrow = 3)

t <- c(2/3, 1)

nSims <- 25000 # 100 # 1000 # 10 #

# Calling the simulation function
out3 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 600,
                   TestStatBin = "UnPooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Bonf", info_frac = t,
                   typeOfDesign = "noEarlyEfficacy",
                   MultipleWinners = F, Selection = T, SelectionScale = "delta",
                   SelectionCriterion = "best", SelectionParmeter = 1,
                   ImplicitSSR = "All",
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out3

# Results:
# AdaptGMCP output:
# FWER = 0.0208, CI = (0.01903, 0.02257)
# East output:
# FWER = 0.02125
#

# Example provided by Cyrus:
# For the binomial example with four treatment arms versus a common control,
# use the following design inputs. Two stage design. alpha=0.025, one sided.
# N=100/arm with interim analysis at 50/arm and LDOF spending function (i.e.
# spend 0.0015253 of the available 0.025 at stage 1).
# pi_0=0.1. pi_1=pi_2=pi_3=pi_4=0.22.
# In East MAMS this has 0.815 power. You can start with this and then add the
# weights and then add the secondary endpoints. For the secondary endpoints,
# use the same values of pi. Assume a correlation of 0.5 between the primary
# and secondary endpoints.

# TEST 2 #########################
# Setting: 1 endpoint and 2 stages
#          Test type: Sidak
#          Early stopping with OF efficacy boundary
#          Multiple winners not allowed (East only has this option)

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = rep(0.22, 5)) # global null hypothesis
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- rep(1/4, 4)
trans_mat <- matrix(c(0, 1/3, 1/3, 1/3,
                      1/3, 0, 1/3, 1/3,
                      1/3, 1/3, 0, 1/3,
                      1/3, 1/3, 1/3, 0), byrow = T, nrow = 4)

t <- c(0.5, 1)

nSims <- 25000 # 100 # 1000 # 10 #

# Calling the simulation function
out2 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                   TestStatBin = "Pooled", # UnPooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Sidak", info_frac = t,
                   typeOfDesign = "asOF",
                   MultipleWinners = F, Selection = T, SelectionScale = "delta",
                   SelectionCriterion = "best", SelectionParmeter = 1,
                   ImplicitSSR = "All",
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out2

# Results:
# AdaptGMCP output:
# FWER = 0.0172, CI = (0.01559, 0.01881)
# East output:
# FWER = 0.0188

# TEST 1 #########################
# Comparing FWER with East
# Setting: 1 endpoint and 2 stages
#          Test type: Sidak
#          No early stopping

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = rep(0.22, 5)) # global null hypothesis
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- rep(1/4, 4)
trans_mat <- matrix(c(0, 1/3, 1/3, 1/3,
                      1/3, 0, 1/3, 1/3,
                      1/3, 1/3, 0, 1/3,
                      1/3, 1/3, 1/3, 0), byrow = T, nrow = 4)

t <- c(0.5, 1)

nSims <- 25000 # 100 # 1000 # 10 #

# Calling the simulation function
out1 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                   TestStatBin = "Pooled", # UnPooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Sidak", info_frac = t,
                   typeOfDesign = "noEarlyEfficacy",
                   MultipleWinners = F, Selection = T,
                   SelectionCriterion = "best", SelectionParmeter = 1,
                   ImplicitSSR = "All",
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out1

# Results:
# AdaptGMCP output:
# FWER = 0.01828, CI = (0.01662, 0.01994)
# East output:
# FWER = 0.02
# East value is just outside the 95% CI computed by AdaptGMCP.


