# File: FS_GMCP_Sim_Bin_Example.R
# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

# Script for testing simMAMSMEP() for fixed sample designs witth a binomial endpoint

library(AdaptGMCP)
library(graphicalMCP)
library(gMCP)
# source("internalData/try_gmcp.R")

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

# TEST 7.2 ###########################
# Comparing FWER with East
# Setting: 1 endpoint and 1 stage (fixed sample design)
#          Test type: Weighted Dunnett

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = rep(0.22, 5))
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- c(0.4, 0.4, 0.1, 0.1)
# To compare with East, edge weights must be set to 0.
trans_mat <- matrix(rep(0, 16), byrow = T, nrow = 4)

# epCorr <- matrix(rep(0.5, 16), byrow = T, nrow = 4)
# diag(epCorr) <- 1

nSims <- 10 # 25000 # 10000 # 100 # 1000 #

# Calling the simulation function
out7 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                   TestStatBin = "UnPooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Dunnett", info_frac = 1,
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = F) # T)

out7

# #########################################################################
# ############################ TESTS FOR POWER ############################
# ############################ FOUR ENDPOINTS #############################
# #########################################################################
#
# # TEST4-4EP-POW ################################
# # Comparing power with East
# # Setting: 4 endpoints and 1 stage (fixed sample design)
# #          Parallel gatekeeping, Test type: Bonferroni
#
# num_arms <- 2 # one treatment arm versus one control
# ep_type <- list(EP1 = "Binary", EP2 = "Binary",
#                 EP3 = "Binary", EP4 = "Binary") # Four binary endpoints
# arm_props <- list(EP1 = c(0.1, 0.22), EP2 = c(0.1, 0.22),
#                   EP3 = c(0.1, 0.22), EP4 = c(0.1, 0.22))
# alloc_ratio <- c(1, 1) # balanced design
#
# # Correlation between endpoints:
# # All 4 endpoints correlated with each other with each pairwise correlation = 0.5
# corr <- matrix(c(1, 0.5, 0.5, 0.5,
#                  0.5, 1, 0.5, 0.5,
#                  0.5, 0.5, 1, 0.5,
#                  0.5, 0.5, 0.5, 1), byrow = T, nrow = 4)
#
# # Graph for parallel gatekeeping:
# # First family: EP1 and EP2
# # Second family: EP3 and EP4
# # Test second family only if at least one hypothesis from the first family
# # is rejected.
# weights <- c(1/2, 1/2, 0, 0)
# trans_mat <- matrix(c(0, 0, 1/2, 1/2,
#                       0, 0, 1/2, 1/2,
#                       0, 0, 0, 1,
#                       0, 0, 1, 0), byrow = T, nrow = 4)
#
# nSims <- 25000 # 10000 # 100 # 1000 # 10 #
#
# # Calling the simulation function
# out4.4EP <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 200,
#                        TestStatBin = "Pooled", # UnPooled",
#                        FWERControl = "CombinationTest",
#                        nArms = num_arms, nEps = 4, lEpType = ep_type,
#                        Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
#                        EP.Corr = corr, WI = weights, G = trans_mat,
#                        test.type = "Bonf", info_frac = 1,
#                        nSimulation = nSims, Seed="Random", EastSumStat = NULL,
#                        Parallel = T)
#
# out4.4EP
#
# # TEST3-4EP-POW ################################
# # Comparing power with East
# # Setting: 4 endpoints and 1 stage (fixed sample design)
# #          Parallel gatekeeping, Test type: Bonferroni
#
# num_arms <- 2 # one treatment arm versus one control
# ep_type <- list(EP1 = "Binary", EP2 = "Binary",
#                 EP3 = "Binary", EP4 = "Binary") # Four binary endpoints
# arm_props <- list(EP1 = c(0.1, 0.22), EP2 = c(0.1, 0.22),
#                   EP3 = c(0.1, 0.1), EP4 = c(0.1, 0.22))
# alloc_ratio <- c(1, 1) # balanced design
#
# # Correlation between endpoints:
# # All 4 endpoints correlated with each other with each pairwise correlation = 0.5
# corr <- matrix(c(1, 0.5, 0.5, 0.5,
#                  0.5, 1, 0.5, 0.5,
#                  0.5, 0.5, 1, 0.5,
#                  0.5, 0.5, 0.5, 1), byrow = T, nrow = 4)
#
# # Graph for parallel gatekeeping:
# # First family: EP1 and EP2
# # Second family: EP3 and EP4
# # Test second family only if at least one hypothesis from the first family
# # is rejected.
# weights <- c(1/2, 1/2, 0, 0)
# trans_mat <- matrix(c(0, 0, 1/2, 1/2,
#                       0, 0, 1/2, 1/2,
#                       0, 0, 0, 1,
#                       0, 0, 1, 0), byrow = T, nrow = 4)
#
# nSims <- 25000 # 10000 # 100 # 1000 # 10 #
#
# # Calling the simulation function
# out3.4EP <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 200,
#                        TestStatBin = "Pooled", # UnPooled",
#                        FWERControl = "CombinationTest",
#                        nArms = num_arms, nEps = 4, lEpType = ep_type,
#                        Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
#                        EP.Corr = corr, WI = weights, G = trans_mat,
#                        test.type = "Bonf", info_frac = 1,
#                        nSimulation = nSims, Seed="Random", EastSumStat = NULL,
#                        Parallel = T)
#
# out3.4EP
#
# # TEST2-4EP-POW ################################
# # Comparing power with East
# # Setting: 4 endpoints and 1 stage (fixed sample design)
# #          Parallel gatekeeping, Test type: Bonferroni
#
# num_arms <- 2 # one treatment arm versus one control
# ep_type <- list(EP1 = "Binary", EP2 = "Binary",
#                 EP3 = "Binary", EP4 = "Binary") # Four binary endpoints
# arm_props <- list(EP1 = c(0.1, 0.22), EP2 = c(0.1, 0.22),
#                   EP3 = c(0.1, 0.1), EP4 = c(0.1, 0.1))
# alloc_ratio <- c(1, 1) # balanced design
#
# # Correlation between endpoints:
# # EP1 is correlated with EP3 and EP2 with EP4. No other correlations.
# corr <- matrix(c(1, 0, 0.5, 0,
#                  0, 1, 0, 0.5,
#                  0.5, 0, 1, 0,
#                  0, 0.5, 0, 1), byrow = T, nrow = 4)
#
# # Graph for parallel gatekeeping:
# # First family: EP1 and EP2
# # Second family: EP3 and EP4
# # Test second family only if at least one hypothesis from the first family
# # is rejected.
# weights <- c(1/2, 1/2, 0, 0)
# trans_mat <- matrix(c(0, 0, 1/2, 1/2,
#                       0, 0, 1/2, 1/2,
#                       0, 0, 0, 1,
#                       0, 0, 1, 0), byrow = T, nrow = 4)
#
# nSims <- 25000 # 10000 # 100 # 1000 # 10 #
#
# # Calling the simulation function
# out2.4EP <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 200,
#                        TestStatBin = "Pooled", # UnPooled",
#                        FWERControl = "CombinationTest",
#                        nArms = num_arms, nEps = 4, lEpType = ep_type,
#                        Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
#                        EP.Corr = corr, WI = weights, G = trans_mat,
#                        test.type = "Bonf", info_frac = 1,
#                        nSimulation = nSims, Seed="Random", EastSumStat = NULL,
#                        Parallel = T)
#
# out2.4EP
#
# # TEST1-4EP-POW ################################
# # Comparing power with East
# # Setting: 4 endpoints and 1 stage (fixed sample design)
# #          Parallel gatekeeping, Test type: Bonferroni
#
# num_arms <- 2 # one treatment arm versus one control
# ep_type <- list(EP1 = "Binary", EP2 = "Binary",
#                 EP3 = "Binary", EP4 = "Binary") # Four binary endpoints
# arm_props <- list(EP1 = c(0.1, 0.22), EP2 = c(0.1, 0.22),
#                   EP3 = c(0.1, 0.1), EP4 = c(0.1, 0.22))
# alloc_ratio <- c(1, 1) # balanced design
#
# # Correlation between endpoints:
# # EP1 is correlated with EP3 and EP2 with EP4. No other correlations.
# corr <- matrix(c(1, 0, 0.5, 0,
#                  0, 1, 0, 0.5,
#                  0.5, 0, 1, 0,
#                  0, 0.5, 0, 1), byrow = T, nrow = 4)
#
# # Graph for parallel gatekeeping:
# # First family: EP1 and EP2
# # Second family: EP3 and EP4
# # Test second family only if at least one hypothesis from the first family
# # is rejected.
# weights <- c(1/2, 1/2, 0, 0)
# trans_mat <- matrix(c(0, 0, 1/2, 1/2,
#                       0, 0, 1/2, 1/2,
#                       0, 0, 0, 1,
#                       0, 0, 1, 0), byrow = T, nrow = 4)
#
# nSims <- 25000 # 10000 # 100 # 1000 # 10 #
#
# # Calling the simulation function
# out1.4EP <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 200,
#                        TestStatBin = "Pooled", # UnPooled",
#                        FWERControl = "CombinationTest",
#                        nArms = num_arms, nEps = 4, lEpType = ep_type,
#                        Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
#                        EP.Corr = corr, WI = weights, G = trans_mat,
#                        test.type = "Bonf", info_frac = 1,
#                        nSimulation = nSims, Seed="Random", EastSumStat = NULL,
#                        Parallel = T)
#
# out1.4EP

#########################################################################
############################ TESTS FOR FWER #############################
############################ FOUR ENDPOINTS #############################
#########################################################################
#
# # TEST1-4EP ################################
# # Comparing FWER with East
# # Setting: 4 endpoints and 1 stage (fixed sample design)
# #          Parallel gatekeeping, Test type: Bonferroni
#
# num_arms <- 2 # one treatment arm versus one control
# ep_type <- list(EP1 = "Binary", EP2 = "Binary",
#                 EP3 = "Binary", EP4 = "Binary") # Four binary endpoints
# arm_props <- list(EP1 = c(0.22, 0.22), EP2 = c(0.22, 0.22),
#                   EP3 = c(0.22, 0.22), EP4 = c(0.22, 0.22))
# alloc_ratio <- c(1, 1) # balanced design
#
# # Correlation between endpoints:
# # EP1 is correlated with EP3 and EP2 with EP4. No other correlations.
# corr <- matrix(c(1, 0, 0.5, 0,
#                  0, 1, 0, 0.5,
#                  0.5, 0, 1, 0,
#                  0, 0.5, 0, 1), byrow = T, nrow = 4)
#
# # Graph for parallel gatekeeping:
# # First family: EP1 and EP2
# # Second family: EP3 and EP4
# # Test second family only if at least one hypothesis from the first family
# # is rejected.
# weights <- c(1/2, 1/2, 0, 0)
# trans_mat <- matrix(c(0, 0, 1/2, 1/2,
#                       0, 0, 1/2, 1/2,
#                       0, 0, 0, 1,
#                       0, 0, 1, 0), byrow = T, nrow = 4)
#
# nSims <- 25000 # 10000 # 100 # 1000 # 10 #
#
# # Calling the simulation function
# out1.4EP <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 200,
#                        TestStatBin = "Pooled", # UnPooled",
#                        FWERControl = "CombinationTest",
#                        nArms = num_arms, nEps = 4, lEpType = ep_type,
#                        Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
#                        EP.Corr = corr, WI = weights, G = trans_mat,
#                        test.type = "Bonf", info_frac = 1,
#                        nSimulation = nSims, Seed="Random", EastSumStat = NULL,
#                        Parallel = T)
#
# out1.4EP

#########################################################################
############################ TESTS FOR POWER ############################
############################ TWO ENDPOINTS ##############################
#########################################################################

# TEST 5-2EP ###########################
# Comparing power with East
# Setting: 2 endpoints and 1 stage (fixed sample design)
#          Test type: Bonferroni

num_arms <- 2 # one treatment arm versus one control
ep_type <- list(EP1 = "Binary", EP2 = "Binary") # Two binary endpoints
arm_props <- list(EP1 = c(0.1, 0.22), EP2 = c(0.1, 0.22))
alloc_ratio <- c(1, 1) # balanced design
corr <- matrix(c(1, 0.5, 0.5, 1), byrow = T, nrow = 2)

weights <- c(1/2, 1/2)
# To compare with East, all edge weights must be set to 0.
trans_mat <- matrix(c(0, 0,
                      0, 0), byrow = T, nrow = 2)
# trans_mat <- matrix(c(0, 1,
#                       1, 0), byrow = T, nrow = 2)

nSims <- 25000 # 10000 # 100 # 1000 # 10 #

# Calling the simulation function
out5.2EP <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 400,
                       TestStatBin = "Pooled", # UnPooled",
                       FWERControl = "CombinationTest",
                       nArms = num_arms, nEps = 2, lEpType = ep_type,
                       Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                       EP.Corr = corr, WI = weights, G = trans_mat,
                       test.type = "Sidak", info_frac = 1,
                       nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                       Parallel = T)

out5.2EP

# TEST 4-2EP ###########################
# Comparing power with East
# Setting: 2 endpoints and 1 stage (fixed sample design)
#          Test type: Bonferroni

num_arms <- 2 # one treatment arm versus one control
ep_type <- list(EP1 = "Binary", EP2 = "Binary") # Two binary endpoints
arm_props <- list(EP1 = c(0.1, 0.22), EP2 = c(0.1, 0.22))
alloc_ratio <- c(1, 1) # balanced design
corr <- matrix(c(1, 0.5, 0.5, 1), byrow = T, nrow = 2)

weights <- c(1/2, 1/2)
# To compare with East, all edge weights must be set to 0.
trans_mat <- matrix(c(0, 0,
                      0, 0), byrow = T, nrow = 2)
# trans_mat <- matrix(c(0, 1,
#                       1, 0), byrow = T, nrow = 2)

nSims <- 25000 # 10000 # 100 # 1000 # 10 #

# Calling the simulation function
out4.2EP <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 400,
                       TestStatBin = "Pooled", # UnPooled",
                       FWERControl = "CombinationTest",
                       nArms = num_arms, nEps = 2, lEpType = ep_type,
                       Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                       EP.Corr = corr, WI = weights, G = trans_mat,
                       test.type = "Bonf", info_frac = 1,
                       nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                       Parallel = T)

out4.2EP

#########################################################################
############################ TESTS FOR FWER #############################
############################ TWO ENDPOINTS ##############################
#########################################################################

# TEST 3-2EP ###########################
# Comparing FWER with East
# Setting: 2 endpoints and 1 stage (fixed sample design)
#          Test type: Sidak

num_arms <- 2 # one treatment arm versus one control
ep_type <- list(EP1 = "Binary", EP2 = "Binary") # Two binary endpoints
arm_props <- list(EP1 = c(0.22, 0.22), EP2 = c(0.22, 0.22))
alloc_ratio <- c(1, 1) # balanced design
corr <- matrix(c(1, 0.5, 0.5, 1), byrow = T, nrow = 2)

weights <- c(1/2, 1/2)
# To compare with East, all edge weights must be set to 0.
trans_mat <- matrix(c(0, 0,
                      0, 0), byrow = T, nrow = 2)
# trans_mat <- matrix(c(0, 1,
#                       1, 0), byrow = T, nrow = 2)

nSims <- 25000 # 10000 # 100 # 1000 # 10 #

# Calling the simulation function
out3.2EP <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 400,
                       TestStatBin = "Pooled", # UnPooled",
                       FWERControl = "CombinationTest",
                       nArms = num_arms, nEps = 2, lEpType = ep_type,
                       Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                       EP.Corr = corr, WI = weights, G = trans_mat,
                       test.type = "Sidak", info_frac = 1,
                       nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                       Parallel = T)

out3.2EP

# TEST 2-2EP ###########################
# Comparing FWER with East
# Setting: 2 endpoints and 1 stage (fixed sample design)
#          Test type: Bonferroni

num_arms <- 2 # one treatment arm versus one control
ep_type <- list(EP1 = "Binary", EP2 = "Binary") # Two binary endpoints
arm_props <- list(EP1 = c(0.1, 0.1), EP2 = c(0.1, 0.1))
alloc_ratio <- c(1, 1) # balanced design
corr <- matrix(c(1, 0.5, 0.5, 1), byrow = T, nrow = 2)

weights <- c(1/2, 1/2)
# To compare with East, all edge weights must be set to 0.
trans_mat <- matrix(c(0, 0,
                      0, 0), byrow = T, nrow = 2)
# trans_mat <- matrix(c(0, 1,
#                       1, 0), byrow = T, nrow = 2)

nSims <- 25000 # 10000 # 100 # 1000 # 10 #

# Calling the simulation function
out2.2EP <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 400,
                       TestStatBin = "Pooled", # UnPooled",
                       FWERControl = "CombinationTest",
                       nArms = num_arms, nEps = 2, lEpType = ep_type,
                       Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                       EP.Corr = corr, WI = weights, G = trans_mat,
                       test.type = "Bonf", info_frac = 1,
                       nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                       Parallel = T)

out2.2EP

# TEST 1-2EP ###########################
# Comparing FWER with East
# Setting: 2 endpoints and 1 stage (fixed sample design)
#          Test type: Bonferroni

num_arms <- 2 # one treatment arm versus one control
ep_type <- list(EP1 = "Binary", EP2 = "Binary") # Two binary endpoints
arm_props <- list(EP1 = c(0.22, 0.22), EP2 = c(0.22, 0.22))
alloc_ratio <- c(1, 1) # balanced design
corr <- matrix(c(1, 0.5, 0.5, 1), byrow = T, nrow = 2)

weights <- c(1/2, 1/2)
# To compare with East, all edge weights must be set to 0.
trans_mat <- matrix(c(0, 0,
                      0, 0), byrow = T, nrow = 2)
# trans_mat <- matrix(c(0, 1,
#                       1, 0), byrow = T, nrow = 2)

nSims <- 25000 # 10000 # 100 # 1000 # 10 #

# Calling the simulation function
out1.2EP <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 400,
                   TestStatBin = "Pooled", # UnPooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 2, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   EP.Corr = corr, WI = weights, G = trans_mat,
                   test.type = "Bonf", info_frac = 1,
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out1.2EP


#########################################################################
############################ SINGLE ENDPOINT ############################
#########################################################################

# TEST 14 ###########################
# Comparing power of the specified scenario with East
# Setting: 1 endpoint and 1 stage (fixed sample design)
#          Test type: Dunett single step

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = c(0.1, 0.22, 0.22, 0.22, 0.22))
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- rep(1/4, 4)
trans_mat <- matrix(rep(0, 16), byrow = T, nrow = 4)

nSims <- 25000 # 10000 # 1000 # 10 # 100 #

# Calling the simulation function
out14 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                    TestStatBin = "Pooled", # UnPooled",
                    FWERControl = "CombinationTest",
                    nArms = num_arms, nEps = 1, lEpType = ep_type,
                    Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                    WI = weights, G = trans_mat,
                    test.type = "Dunnett", info_frac = 1,
                    nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                    Parallel = T)

out14

# TEST 13 ###########################
# Comparing power of the specified scenario with East
# Setting: 1 endpoint and 1 stage (fixed sample design)
#          Test type: Fallback (Bonferroni with special weights in AdaptGMCP)
# We are using Bonferroni here only to trick AdaptGMCP into doing this test.

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = c(0.1, 0.22, 0.22, 0.22, 0.22))
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- c(0.4, 0.3, 0.1, 0.1)
# Note that the node weights add up to less than 1, which is allowed.

trans_mat <- matrix(c(0, 1, 0, 0,
                      0, 0, 1, 0,
                      0, 0, 0, 1,
                      0, 0, 0, 0), byrow = T, nrow = 4)

nSims <- 25000 # 10000 # 1000 # 10 # 100 #

# Calling the simulation function
out13 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                    TestStatBin = "Pooled", # UnPooled",
                    FWERControl = "CombinationTest",
                    nArms = num_arms, nEps = 1, lEpType = ep_type,
                    Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                    WI = weights, G = trans_mat,
                    test.type = "Bonf", info_frac = 1,
                    nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                    Parallel = T)

out13

########## Testing this same example using graphicalMCP
weights <- c(0.4, 0.3, 0.1, 0.1)
hyp_names <- c("H1", "H2", "H3", "H4")
gr <- fallback(weights)
plot(gr)

arm_props <- list(EP1 = c(0.1, 0.22, 0.22, 0.22, 0.22))
alpha <- 0.025
ss <- 500
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

p <- unlist(arm_props)
unpooled_var <-
  p[-1] * (1-p[-1]) / (ss/5) + p[1] * (1-p[1]) / (ss/5)
non_centr_par <- (p[-1] - p[1]) / sqrt(unpooled_var)

marg_pow <- pnorm(
  qnorm(alpha, lower.tail = F),
  mean = non_centr_par,
  sd = 1,
  lower.tail = F)
names(marg_pow) <- hyp_names
marg_pow

corr <- matrix(c(1, 0.5, 0.5, 0.5,
                 0.5, 1, 0.5, 0.5,
                 0.5, 0.5, 1, 0.5,
                 0.5, 0.5, 0.5, 1), byrow = T, nrow = 4)
colnames(corr) <- rownames(corr) <- hyp_names
corr

pow13 <- graph_calculate_power(gr, alpha = alpha, sim_corr = corr,
                             sim_n = 1e5, power_marginal = marg_pow)

pow13

######### Testing the example using gMCP


# Another run with node weights adding up to 1.
weights <- c(0.4, 0.3, 0.2, 0.1)

out13.2 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                    TestStatBin = "Pooled", # UnPooled",
                    FWERControl = "CombinationTest",
                    nArms = num_arms, nEps = 1, lEpType = ep_type,
                    Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                    WI = weights, G = trans_mat,
                    test.type = "Bonf", info_frac = 1,
                    nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                    Parallel = T)

out13.2

# Testing this same example using graphicalMCP
weights <- c(0.4, 0.3, 0.2, 0.1)
gr <- fallback(weights)
plot(gr)

pow13.2 <- graph_calculate_power(gr, alpha = alpha, sim_corr = corr,
                               sim_n = 1e5, power_marginal = marg_pow)

pow13.2

# TEST 12 ###########################
# Comparing power of the specified scenario with East
# Setting: 1 endpoint and 1 stage (fixed sample design)
#          Test type: Fixed Sequence (Bonferroni with special weights in AdaptGMCP)
# We are using Bonferroni here only to trick AdaptGMCP into doing a fixed
# sequence test.

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = c(0.1, 0.22, 0.22, 0.22, 0.22))
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- c(1, 0, 0, 0)
trans_mat <- matrix(c(0, 1, 0, 0,
                      0, 0, 1, 0,
                      0, 0, 0, 1,
                      0, 0, 0, 0), byrow = T, nrow = 4)

nSims <- 25000 # 10000 # 1000 # 10 # 100 #

# Calling the simulation function
out12 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                    TestStatBin = "Pooled", # UnPooled",
                    FWERControl = "CombinationTest",
                    nArms = num_arms, nEps = 1, lEpType = ep_type,
                    Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                    WI = weights, G = trans_mat,
                    test.type = "Bonf", info_frac = 1,
                    nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                    Parallel = T)

out12

# Testing this same example using graphicalMCP
hyp_names <- c("H1", "H2", "H3", "H4")
gr <- fixed_sequence(4)
plot(gr)

arm_props <- list(EP1 = c(0.1, 0.22, 0.22, 0.22, 0.22))
alpha <- 0.025
ss <- 500
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

p <- unlist(arm_props)
unpooled_var <-
  p[-1] * (1-p[-1]) / (ss/5) + p[1] * (1-p[1]) / (ss/5)
non_centr_par <- (p[-1] - p[1]) / sqrt(unpooled_var)

marg_pow <- pnorm(
  qnorm(alpha, lower.tail = F),
  mean = non_centr_par,
  sd = 1,
  lower.tail = F)
names(marg_pow) <- hyp_names
marg_pow

corr <- matrix(c(1, 0.5, 0.5, 0.5,
                 0.5, 1, 0.5, 0.5,
                 0.5, 0.5, 1, 0.5,
                 0.5, 0.5, 0.5, 1), byrow = T, nrow = 4)
colnames(corr) <- rownames(corr) <- hyp_names
corr

pow <- graph_calculate_power(gr, alpha = alpha, sim_corr = corr,
                             sim_n = 1e5, power_marginal = marg_pow)

pow


# TEST 11 ###########################
# Comparing power of the specified scenario with East
# Setting: 1 endpoint and 1 stage (fixed sample design)
#          Test type: Bonferroni (Holm in East)

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = c(0.1, 0.22, 0.22, 0.22, 0.22))
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- rep(1/4, 4)
# To compare with East, all edge weights must be set to 1/3.
trans_mat <- matrix(c(0, 1/3, 1/3, 1/3,
                      1/3, 0, 1/3, 1/3,
                      1/3, 1/3, 0, 1/3,
                      1/3, 1/3, 1/3, 0), byrow = T, nrow = 4)

nSims <- 25000 # 10000 # 1000 # 10 # 100 #

# Calling the simulation function
out11 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                   TestStatBin = "Pooled", # UnPooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Bonf", info_frac = 1,
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out11

# Testing this same example using graphicalMCP
hyp_names <- c("H1", "H2", "H3", "H4")
gr <- bonferroni_holm(4)
plot(gr)

arm_props <- list(EP1 = c(0.1, 0.22, 0.22, 0.22, 0.22))
alpha <- 0.025
ss <- 500
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

p <- unlist(arm_props)
unpooled_var <-
  p[-1] * (1-p[-1]) / (ss/5) + p[1] * (1-p[1]) / (ss/5)
non_centr_par <- (p[-1] - p[1]) / sqrt(unpooled_var)

marg_pow <- pnorm(
  qnorm(alpha, lower.tail = F),
  mean = non_centr_par,
  sd = 1,
  lower.tail = F)
names(marg_pow) <- hyp_names
marg_pow

corr <- matrix(c(1, 0.5, 0.5, 0.5,
                 0.5, 1, 0.5, 0.5,
                 0.5, 0.5, 1, 0.5,
                 0.5, 0.5, 0.5, 1), byrow = T, nrow = 4)
colnames(corr) <- rownames(corr) <- hyp_names
corr

pow11 <- graph_calculate_power(gr, alpha = alpha, sim_corr = corr,
                             sim_n = 1e5, power_marginal = marg_pow)

pow11

# TEST 10 ###########################
# Comparing power of the specified scenario with East
# Setting: 1 endpoint and 1 stage (fixed sample design)
#          Test type: Weighted Bonferroni

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = c(0.1, 0.22, 0.22, 0.22, 0.22))
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- c(0.4, 0.4, 0.1, 0.1)
# To compare with East, all edge weights must be set to 0.
trans_mat <- matrix(rep(0, 16), byrow = T, nrow = 4)
# trans_mat <- matrix(c(0, 1/3, 1/3, 1/3,
#                       1/3, 0, 1/3, 1/3,
#                       1/3, 1/3, 0, 1/3,
#                       1/3, 1/3, 1/3, 0), byrow = T, nrow = 4)

nSims <- 25000 # 10000 # 10 # 100 # 1000 #

# Calling the simulation function
out10 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                   TestStatBin = "Pooled", # UnPooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Bonf", info_frac = 1,
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out10

# Testing using graphicalMCP:
hyp_names <- c("H1", "H2", "H3", "H4")
weights <- c(0.4, 0.4, 0.1, 0.1)
gr <- bonferroni_weighted(weights)
plot(gr)

arm_props <- list(EP1 = c(0.1, 0.22, 0.22, 0.22, 0.22))
alpha <- 0.025
ss <- 500
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

p <- unlist(arm_props)
unpooled_var <-
  p[-1] * (1-p[-1]) / (ss/5) + p[1] * (1-p[1]) / (ss/5)
non_centr_par <- (p[-1] - p[1]) / sqrt(unpooled_var)

marg_pow <- pnorm(
  qnorm(alpha, lower.tail = F),
  mean = non_centr_par,
  sd = 1,
  lower.tail = F)
names(marg_pow) <- hyp_names

marg_pow

corr <- matrix(c(1, 0.5, 0.5, 0.5,
                 0.5, 1, 0.5, 0.5,
                 0.5, 0.5, 1, 0.5,
                 0.5, 0.5, 0.5, 1), byrow = T, nrow = 4)
colnames(corr) <- rownames(corr) <- hyp_names
corr

pow <- graph_calculate_power(gr, alpha = alpha, sim_corr = corr,
                             sim_n = 1e5, power_marginal = marg_pow)

pow

# TEST 9 ###########################
# Comparing power of the specified scenario with East
# Setting: 1 endpoint and 1 stage (fixed sample design)
#          Test type: Sidak

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = c(0.1, 0.22, 0.22, 0.22, 0.22))
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- rep(1/4, 4)
# To compare with East, all edge weights must be set to 0.
trans_mat <- matrix(rep(0, 16), byrow = T, nrow = 4)
# trans_mat <- matrix(c(0, 1/3, 1/3, 1/3,
#                       1/3, 0, 1/3, 1/3,
#                       1/3, 1/3, 0, 1/3,
#                       1/3, 1/3, 1/3, 0), byrow = T, nrow = 4)

nSims <- 25000 # 10000 # 10 # 1000 # 100 #

# Calling the simulation function
out9 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                  TestStatBin = "Pooled", # UnPooled",
                  FWERControl = "CombinationTest",
                  nArms = num_arms, nEps = 1, lEpType = ep_type,
                  Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                  WI = weights, G = trans_mat,
                  test.type = "Sidak", info_frac = 1,
                  nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                  Parallel = T)

out9

# Testing this same example using graphicalMCP
hyp_names <- c("H1", "H2", "H3", "H4")
gr <- sidak(4)
plot(gr)

arm_props <- list(EP1 = c(0.1, 0.22, 0.22, 0.22, 0.22))
alpha <- 0.025
ss <- 500
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

p <- unlist(arm_props)
unpooled_var <-
  p[-1] * (1-p[-1]) / (ss/5) + p[1] * (1-p[1]) / (ss/5)
non_centr_par <- (p[-1] - p[1]) / sqrt(unpooled_var)

marg_pow <- pnorm(
  qnorm(alpha, lower.tail = F),
  mean = non_centr_par,
  sd = 1,
  lower.tail = F)
names(marg_pow) <- hyp_names
marg_pow

corr <- matrix(c(1, 0.5, 0.5, 0.5,
                 0.5, 1, 0.5, 0.5,
                 0.5, 0.5, 1, 0.5,
                 0.5, 0.5, 0.5, 1), byrow = T, nrow = 4)
colnames(corr) <- rownames(corr) <- hyp_names
corr

pow9 <- graph_calculate_power(gr, alpha = alpha, sim_corr = corr,
                               sim_n = 1e5, power_marginal = marg_pow)

pow9

# TEST 8 ###########################
# Comparing power of the specified scenario with East
# Setting: 1 endpoint and 1 stage (fixed sample design)
#          Test type: Bonferroni

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = c(0.1, 0.22, 0.22, 0.22, 0.22))
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- rep(1/4, 4)
# To compare with East, all edge weights must be set to 0.
trans_mat <- matrix(rep(0, 16), byrow = T, nrow = 4)
# trans_mat <- matrix(c(0, 1/3, 1/3, 1/3,
#                       1/3, 0, 1/3, 1/3,
#                       1/3, 1/3, 0, 1/3,
#                       1/3, 1/3, 1/3, 0), byrow = T, nrow = 4)

nSims <- 100 # 25000 # 10000 # 1000 # 10 #

# Calling the simulation function
out8 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                  TestStatBin = "Pooled", # UnPooled",
                  FWERControl = "CombinationTest",
                  nArms = num_arms, nEps = 1, lEpType = ep_type,
                  Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                  WI = weights, G = trans_mat,
                  test.type = "Bonf", info_frac = 1,
                  nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                  Parallel = T)

out8

######### Testing the example using gMCP
weights <- rep(1/4, 4)
trans_mat <- matrix(rep(0, 16), byrow = T, nrow = 4)
hyp_names <- c("H1", "H2", "H3", "H4")

# Creating a new gMCP graph using these weights and transition matrix
gGr <- new("graphMCP", m=trans_mat, weights=weights)

p <- c(0.1, 0.22, 0.22, 0.22, 0.22)
n <- c(100, 100, 100, 100, 100)
se <- sqrt(p[-1]*(1-p[-1])/n[-1] + p[1]*(1-p[1])/n[1])

# Calculating the z-scores / non-centrality parameters to be passed as means
z_means <- (p[-1]-p[1])/se

# Correlation matrix to be used for MVN data generation
mCr <- matrix(c(1, 0.5, 0.5, 0.5,
                0.5, 1, 0.5, 0.5,
                0.5, 0.5, 1, 0.5,
                0.5, 0.5, 0.5, 1), byrow = T, nrow = 4)

# Calling gMCP
pow8 <- gMCP::calcPower(mean = z_means, corr.sim = mCr,
                        graph = gGr, alpha = 0.025, n.sim = 10000)
pow8

############ Testing this same example using graphicalMCP
weights <- rep(1/4, 4)
trans_mat <- matrix(rep(0, 16), byrow = T, nrow = 4)
hyp_names <- c("H1", "H2", "H3", "H4")
gr <- graph_create(weights, trans_mat, hyp_names)
plot(gr)

# The package also provides a shortcut to create the Bonferroni graph:
# gr <- bonferroni(4)

arm_props <- list(EP1 = c(0.1, 0.22, 0.22, 0.22, 0.22))
alpha <- 0.025
ss <- 500
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

p <- unlist(arm_props)
nc <- nt <- ss/5
pooled_pi <- (nt * p[-1] + nc * p[1]) / (nt + nc)
pooled_var <- pooled_pi * (1 - pooled_pi) * (1/nt + 1/nc)
#
unpooled_var <-
  p[-1] * (1-p[-1]) / (ss/5) + p[1] * (1-p[1]) / (ss/5)
non_centr_par <- (p[-1] - p[1]) / sqrt(unpooled_var)
# non_centr_par <- (p[-1] - p[1]) / sqrt(pooled_var)

marg_pow <- pnorm(
  qnorm(alpha, lower.tail = F),
  mean = non_centr_par,
  sd = 1,
  lower.tail = F)
names(marg_pow) <- hyp_names

corr <- matrix(c(1, 0.5, 0.5, 0.5,
                 0.5, 1, 0.5, 0.5,
                 0.5, 0.5, 1, 0.5,
                 0.5, 0.5, 0.5, 1), byrow = T, nrow = 4)
colnames(corr) <- rownames(corr) <- hyp_names

pow <- graph_calculate_power(gr, alpha = alpha, sim_corr = corr,
                             sim_n = 1e5, power_marginal = marg_pow)

pow


# TEST 7 ###########################
# Comparing FWER with East
# Setting: 1 endpoint and 1 stage (fixed sample design)
#          Test type: Weighted Sidak

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = rep(0.22, 5))
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- c(0.4, 0.4, 0.1, 0.1)
# To compare with East, edge weights must be set to 0.
trans_mat <- matrix(rep(0, 16), byrow = T, nrow = 4)
# trans_mat <- matrix(c(0, 1/3, 1/3, 1/3,
#                       1/3, 0, 1/3, 1/3,
#                       1/3, 1/3, 0, 1/3,
#                       1/3, 1/3, 1/3, 0), byrow = T, nrow = 4)

nSims <- 25000 # 10000 # 10 # 100 # 1000 #

# Calling the simulation function
out7 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                   TestStatBin = "Pooled", # UnPooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Sidak", info_frac = 1,
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out7

# TEST 6 ###########################
# Comparing FWER with East
# Setting: 1 endpoint and 1 stage (fixed sample design)
#          Test type: Weighted Bonferroni

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = rep(0.22, 5))
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- c(0.4, 0.4, 0.1, 0.1)
# To compare with East, edge weights must be set to 0.
trans_mat <- matrix(rep(0, 16), byrow = T, nrow = 4)
# trans_mat <- matrix(c(0, 1/3, 1/3, 1/3,
#                       1/3, 0, 1/3, 1/3,
#                       1/3, 1/3, 0, 1/3,
#                       1/3, 1/3, 1/3, 0), byrow = T, nrow = 4)

nSims <- 25000 # 10000 # 10 # 100 # 1000 #

# Calling the simulation function
out6 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                   TestStatBin = "Pooled", # UnPooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Bonf", info_frac = 1,
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out6

######### Testing the example using gMCP
weights <- c(0.4, 0.4, 0.1, 0.1)
trans_mat <- matrix(rep(0, 16), byrow = T, nrow = 4)

# Creating a new gMCP graph using these weights and transition matrix
gGr <- new("graphMCP", m=trans_mat, weights=weights)

p <- rep(0.22, 5)
n <- c(100, 100, 100, 100, 100)
se <- sqrt(p[-1]*(1-p[-1])/n[-1] + p[1]*(1-p[1])/n[1])

# Calculating the z-scores / non-centrality parameters to be passed as means
z_means <- (p[-1]-p[1])/se

# Correlation matrix to be used for MVN data generation
mCr <- matrix(c(1, 0.5, 0.5, 0.5,
                0.5, 1, 0.5, 0.5,
                0.5, 0.5, 1, 0.5,
                0.5, 0.5, 0.5, 1), byrow = T, nrow = 4)

# Calling gMCP
pow6 <- gMCP::calcPower(mean = z_means, corr.sim = mCr,
                        graph = gGr, alpha = 0.025, n.sim = 100000)
pow6

# TEST 5 ###########################
# Comparing FWER with East
# Setting: 1 endpoint and 1 stage (fixed sample design)
#          Test type: Bonferroni

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = rep(0.22, 5))
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- rep(1/4, 4) # c(0.4, 0.4, 0.1, 0.1)
# To compare with East, edge weights must be set to 0.
trans_mat <- matrix(rep(0, 16), byrow = T, nrow = 4)
# trans_mat <- matrix(c(0, 1/3, 1/3, 1/3,
#                       1/3, 0, 1/3, 1/3,
#                       1/3, 1/3, 0, 1/3,
#                       1/3, 1/3, 1/3, 0), byrow = T, nrow = 4)

nSims <- 25000 # 10000 # 10 # 100 # 1000 #

# Calling the simulation function
out5 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                   TestStatBin = "Pooled", # UnPooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Bonf", info_frac = 1,
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out5

# TEST 4 ###########################
# Comparing FWER with East
# Setting: 1 endpoint and 1 stage (fixed sample design)
#          Test type: Weighted Bonferroni

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = c(0.1, 0.1, 0.1, 0.1, 0.1))
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- c(0.4, 0.4, 0.1, 0.1)
# To compare with East, all edge weights must be set to 0.
trans_mat <- matrix(rep(0, 16), byrow = T, nrow = 4)
# trans_mat <- matrix(c(0, 1/3, 1/3, 1/3,
#                       1/3, 0, 1/3, 1/3,
#                       1/3, 1/3, 0, 1/3,
#                       1/3, 1/3, 1/3, 0), byrow = T, nrow = 4)

nSims <- 25000 # 10000 # 10 # 100 # 1000 #

# Calling the simulation function
out4 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                   TestStatBin = "Pooled", # UnPooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Bonf", info_frac = 1,
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out4

# TEST 3 ###########################
# Comparing FWER with East
# Setting: 1 endpoint and 1 stage (fixed sample design)
#          Test type: Simes

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = c(0.1, 0.1, 0.1, 0.1, 0.1))
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- rep(1/4, 4)
trans_mat <- matrix(c(0, 1/3, 1/3, 1/3,
                      1/3, 0, 1/3, 1/3,
                      1/3, 1/3, 0, 1/3,
                      1/3, 1/3, 1/3, 0), byrow = T, nrow = 4)

nSims <- 10000 # 10 # 100 # 1000 #

# Calling the simulation function
out3 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                   TestStatBin = "Pooled", # UnPooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Simes", info_frac = 1,
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out3

# Testing without any alpha propogation
trans_mat <- matrix(rep(0, 16), byrow = T, nrow = 4)
out3.2 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                   TestStatBin = "Pooled", # UnPooled",
                   FWERControl = "CombinationTest",
                   nArms = num_arms, nEps = 1, lEpType = ep_type,
                   Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                   WI = weights, G = trans_mat,
                   test.type = "Simes", info_frac = 1,
                   nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                   Parallel = T)

out3.2$Overall_Powers

# TEST 2 ###########################
# Comparing FWER with East
# Setting: 1 endpoint and 1 stage (fixed sample design)
#          Test type: Sidak

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = c(0.1, 0.1, 0.1, 0.1, 0.1))
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- rep(1/4, 4)
# To compare with East, all edge weights must be set to 0.
trans_mat <- matrix(rep(0, 16), byrow = T, nrow = 4)
# trans_mat <- matrix(c(0, 1/3, 1/3, 1/3,
#                       1/3, 0, 1/3, 1/3,
#                       1/3, 1/3, 0, 1/3,
#                       1/3, 1/3, 1/3, 0), byrow = T, nrow = 4)

nSims <- 10000 # 10 # 100 # 1000 #

# Calling the simulation function
out2 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                  TestStatBin = "Pooled", # UnPooled",
                  FWERControl = "CombinationTest",
                  nArms = num_arms, nEps = 1, lEpType = ep_type,
                  Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                  WI = weights, G = trans_mat,
                  test.type = "Sidak", info_frac = 1,
                  nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                  Parallel = F)

out2

# TEST 1 ###########################
# Comparing FWER with East
# Setting: 1 endpoint and 1 stage (fixed sample design)
#          Test type: Bonferroni

num_arms <- 5 # four treatment arms versus a common control
ep_type <- list(EP1 = "Binary") # One binary endpoint
arm_props <- list(EP1 = c(0.1, 0.1, 0.1, 0.1, 0.1))
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

weights <- rep(1/4, 4)
# To compare with East Bonferroni, all edge weights must be specified as 0.
trans_mat <- matrix(rep(0, 16), byrow = T, nrow = 4)
# trans_mat <- matrix(c(0, 1/3, 1/3, 1/3,
#                       1/3, 0, 1/3, 1/3,
#                       1/3, 1/3, 0, 1/3,
#                       1/3, 1/3, 1/3, 0), byrow = T, nrow = 4)

nSims <- 10000 # 100 # 1000 # 10 #

# Calling the simulation function
out1 <- simMAMSMEP(Method = "CombPValue", alpha = 0.025, SampleSize = 500,
                  TestStatBin = "Pooled", # UnPooled",
                  FWERControl = "CombinationTest",
                  nArms = num_arms, nEps = 1, lEpType = ep_type,
                  Arms.Prop = arm_props, Arms.alloc.ratio = alloc_ratio,
                  WI = weights, G = trans_mat,
                  test.type = "Bonf", info_frac = 1,
                  nSimulation = nSims, Seed="Random", EastSumStat = NULL,
                  Parallel = T)

out1

# Testing this same example using graphicalMCP
weights <- rep(1/4, 4)
trans_mat <- matrix(rep(0, 16), byrow = T, nrow = 4)
hyp_names <- c("H1", "H2", "H3", "H4")
gr <- graph_create(weights, trans_mat, hyp_names)
plot(gr)

arm_props <- list(EP1 = c(0.1, 0.1, 0.1, 0.1, 0.1))
alpha <- 0.025
ss <- 500
alloc_ratio <- c(1, 1, 1, 1, 1) # balanced design

p <- unlist(arm_props)
unpooled_var <-
  p[-1] * (1-p[-1]) / (ss/5) + p[1] * (1-p[1]) / (ss/5)
non_centr_par <- -(p[-1] - p[1]) / sqrt(unpooled_var)

marg_pow <- pnorm(
  qnorm(alpha, lower.tail = F),
  mean = non_centr_par,
  sd = 1,
  lower.tail = F)
names(marg_pow) <- hyp_names

corr <- matrix(c(1, 0.5, 0.5, 0.5,
                 0.5, 1, 0.5, 0.5,
                 0.5, 0.5, 1, 0.5,
                 0.5, 0.5, 0.5, 1), byrow = T, nrow = 4)
colnames(corr) <- rownames(corr) <- hyp_names

pow <- graph_calculate_power(gr, alpha = alpha, sim_corr = corr,
                             sim_n = 1e5, power_marginal = marg_pow)


