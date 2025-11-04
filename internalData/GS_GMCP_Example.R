# GS_GMCP_Example.R
# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

############################################################################
#####   Script to Run Examples for group sequential GMCP Tests  ############
############################################################################

library(AdaptGMCP)

# EXAMPLE 3 #############################################
# Ref: Graph Based MAMS.pdf
# Trial example on page 2
# 2 treatment arms x 2 endpoints, common control: 4 hypo
# Balanced randomization
# Graph as in Fig 1 on page 5
# Interim look at t=0.5
# alpha = 0.025, LDOF spf

# Setting input parameters for the function

# Weights
wi <- c(0.5, 0.5, 0, 0)

# Transition matrix
g <- matrix(c(0, 0.5, 0.5, 0,
              0.5, 0, 0, 0.5,
              0, 1, 0, 0,
              1, 0, 0, 0), byrow = T, nrow = 4)

# Test type
test <- "Dunnett" # "Bonf" # "Partly-Parametric" #

# Type I error
alp <- 0.025

# Info fraction
t <- c(0.5, 1) # This is a fixed sample design

# Design type
des <- "asOF"

# Correlation matrix between test stats
corr <- matrix(c(1, 0.5, 0.5, 0.5,
                 0.5, 1, 0.5, 0.5,
                 0.5, 0.5, 1, 0.5,
                 0.5, 0.5, 0.5, 1), byrow = T, nrow = 4)

# Calling the function
# Using the correlation matrix created above
# Enter these raw p-values when asked:
# Raw p-values for look 1: p1=0.00045, p2=0.0952, p3=0.0225, p4=0.1104
# Look 1 result: H1 rejected, graph updated for H2, H3, H4
# Raw p-values for look 2: p2=0.1121, p3=0.0112, p4=0.1153
# out <- adaptGMCP_PC(WI=wi, G=g, test.type = test, alpha = alp, info_frac = t,
#                     typeOfDesign = des, Correlation = corr, MultipleWinners = T,
#                     Selection = F, UpdateStrategy = F, plotGraphs = T)
# WARNING: The output of this does not match that from the paper very well because
# the output from the paper seems to have been calculated using the default
# correlation matrix rather than what I have created above.

# Checking the output with the default correlation matrix keeping everything
# else the same
out <- adaptGMCP_PC(WI=wi, G=g, test.type = test, alpha = alp, info_frac = t,
                    typeOfDesign = des, MultipleWinners = T,
                    Selection = F, UpdateStrategy = F, plotGraphs = T)
# Now look 1 output matches the paper perfectly.
# WARNING: However, adjusted p-val for H4 at stage 2 does not match that from
# the paper very well (0.1104 by the package vs 0.1153 in the paper).
# Hence the final combined p-val for H4 differs too (package's 0.0416 vs
# paper's 0.0433).

# Now trying out the weight modification after stage 1 as mentioned in the
# example on page 15 of the paper. Everything else the same as before.
out <- adaptGMCP_PC(WI=wi, G=g, test.type = test, alpha = alp, info_frac = t,
                    typeOfDesign = des, MultipleWinners = T,
                    Selection = F, UpdateStrategy = T, plotGraphs = T)
# Output matches with that from tables 6 and 7 on pages 16 and 17 from the paper

# EXAMPLE 4 #############################################
# Ref: Graph Based MAMS.pdf
# Trial example on page 2
# 2 treatment arms x 2 endpoints, common control: 4 hypo
# Balanced randomization
# Total sample size = 210 (70 per arm - control, dose 1, dose 2)
# Graph as in Fig 1 on page 5
# Interim look at t=0.5
# alpha = 0.025, LDOF spf
# Raw p-values and rest of the details starting from the last para of page 24
# up to the first para on page 27 of the paper.

# Setting input parameters for the function

ss <- 210 # Total sample size
arms <- 3 # Control, dose 1, dose 2
endpts <- 2 # primary endpoint, secondary endpoint
endpt_type <- list(EP1 = "Continuous", EP2 = "Continuous")

# Arm-Wise planned std dev for each endpoint
# Not required for test.type = 'Non-Parametric'
sig <- list(EP1 = c(1,1,1), EP2 = c(1,1,1))

# allocation ratios for the 3 arms
r <- c(1, 1, 1)

wi <- c(0.5, 0.5, 0, 0) # Initial weights for the 4 hypo
g <- matrix(c(0, 1/2, 1/2, 0,
              1/2, 0, 0, 1/2,
              0, 1, 0, 0,
              1, 0, 0, 0), byrow = T, nrow = 4)

test <- "Parametric" # MCP test type

alp <- 0.025 # total type I err

t <- c(0.5, 1) # info fractions for the interim and final look

des <- "asOF" # O'Brien Flemming spending function

adapt <- TRUE
graph <- TRUE

# Specify the following input interactively as prompted:
# Raw p-values at look 1: p1=0.00045, p2=0.0952, p3=0.0225, p4=0.1104
# Hypo to be retained for stage 2: H2,H4 (H3 dropped)
# Modified cumulative SS for Control and Treat2: 88, 87 (incrementally 53 and 53)
# Changed testing strategy for stage 2:
#                           new weights for H2, H4: 0.5, 0.5
#                           edge weights for the 2 edges: 1 each
# Incremental raw p-values for stage 2: p2=0.0299, p4=0.0586
out <- adaptGMCP_CER(nArms=arms, nEps=endpts, SampleSize = ss,
                     EpType = endpt_type, sigma = sig, allocRatio = r,
                     WI = wi, G=g, test.type = test, alpha = alp,
                     info_frac = t, typeOfDesign = des,
                     AdaptStage2 = adapt, plotGraphs = graph)
# Result: outputs matches with the paper

# Example from second para on page 27 of the paper:
# Look 1 input as before
# All hypo to be retained including H3 for stage 2
# No sample size change for the arms
# Changed testing strategy for stage 2:
#                           new weights for H2, H3, H4: 0.5, 0.25, 0.25
#                           edge weights for H2, H3, H4 as per this matrix:
#                           matrix(c(0, 1/3, 2/3,
#                                    1, 0, 0,
#                                    1/2, 1/2, 0), byrow=T, nrow=3)
# Incremental raw p-values for stage 2: p2=0.0299, p3=0.0225, p4=0.0586
out <- adaptGMCP_CER(nArms=arms, nEps=endpts, SampleSize = ss,
                     EpType = endpt_type, sigma = sig, allocRatio = r,
                     WI = wi, G=g, test.type = test, alpha = alp,
                     info_frac = t, typeOfDesign = des,
                     AdaptStage2 = adapt, plotGraphs = graph)
# Result: outputs matches with the paper - compare the numbers in the first
# para on page 28 with the Adapt_Boundary for the intersection hypo H234 of the
# table $Adapt_Test_Tables$Adjusted_Boundary.

#######################################################
# Reference paper: "A Comparison of Two Methods for Adaptive Multi-Arm Two-Stage
# Design", May 2025, Mehta & Kappler
# File: Stagewise vs Cumulative MAMS

# EXAMPLE 1 ################################
# From page 9 of above paper
# Two doses of experimental drug compared with one control, one endpoint.
# Trial will be successful if efficacy could be demonstrated with either dose.
# Thus, there are 2 hypotheses H1 and H2, and there is no hierarchy among them.
# Total SS=210. Balanced randomization of total 70 subjects per arm.
# One interim look at t=0.5. Thus, n=35 for each arm at look 1 and look 2.
# Total alpha=0.025
# LD-OF alpha spending function. Thus, alpha_1=0.00153
################ Stagewise MAMS

# LOOK 1 ANALYSIS >>>>>>>>>>>>>>
# Equal weights to H1 and H2
wi <- c(0.5, 0.5) # weights for H1 and H2 in the graph
g <- matrix(c(0, 1, 1, 0), byrow = T, nrow = 2) # transition matrix for graph
tt <- "Dunnett" # MCP test
err_1 <- 0.025 # total type I error
t <- c(0.5, 1) # info fractions for looks 1 and 2
des <- "asOF"
corr <- matrix(c(1, 0.5, 0.5, 1), byrow = T, nrow = 2) # correlation matrix
                                                       # between test stats
multwin <- T # multiple winners option required
sel <- F # no treatment selection at look 1
upd <- F # no modification of weights and testing strategy at look 1
graph <- T # plot intermediate graphs

# Calling function
# Raw p-value for look 1: p11=0.0294, p21=0.0463
adaptGMCP_PC(WI=wi, G=g, test.type = tt, alpha = err_1, info_frac = t,
             typeOfDesign = des, Correlation = corr, MultipleWinners = multwin,
             Selection = sel, UpdateStrategy = upd, plotGraphs = graph)

# Since no hypothesis is rejected, the trial proceeds to look 2.
# LOOK 2 ANALYSIS >>>>>>>>>>>>>>
# Raw p-value for look 2: p1(2)=0.0475, p2(2)=0.1352
# Choose appropriate option after executing the above adaptGMCP_PC() call
# and specify the raw p-values. Keep the original correlation matrix and other
# input unchanged.

# EXAMPLE 2 ############################
# Cumulative MAMS
# From page 10 of the paper

# Setting input parameters
arms <- 3 # Assuming this includes control arm
endpts <- 1 # number of endpoints
SS <- 210 # total sample size
ep_type <- list(EP1="Continuous")
sig <- list(EP1=c(1,1,1))
comm_stddev <- T

# Ignoring param proc.ctr as it is not applicable here

alloc <- c(1, 1, 1)
wi <- c(0.5, 0.5) # weights for H1 and H2 in the graph
g <- matrix(c(0, 1, 1, 0), byrow = T, nrow = 2) # transition matrix for graph
tt <- "Parametric" # For weighted Dunnett's test
err_1 <- 0.025 # total type I error
t <- c(0.5, 1) # info fractions for looks 1 and 2
des <- "asOF"

# Ignoring params deltaWT, deltaPT1, and gammaA as they are not applicable here

adapt <- T
graph <- T # plot intermediate graphs

# Calling function
# Raw p-value for look 1: p11=0.0294, p21=0.0463
adaptGMCP_CER(nArms = arms, nEps = endpts, SampleSize = SS, EpType = ep_type,
              sigma = sig, CommonStdDev = comm_stddev, allocRatio = alloc,
              WI = wi, G=g, test.type = tt, alpha = err_1, info_frac = t,
              typeOfDesign = des, AdaptStage2 = adapt, plotGraphs = graph)

# WARNING: This function prints p-value scale boundaries for stage 1 and 2
# whereas the paper shows Z-scale boundaries. The latter can be calculated by
# inverting the p-value scale boundaries using the function qnorm(). Set lower
# tail to FALSE since this is a right tailed case.
# Z-scale boundaries computed this way match those from the paper for H1 and H2
# very well, but for H12 they match only for the first 2-3 decimal places.
# Similarly the CERs for H1 and H2 match very well with the paper, but that
# for H12 matches for first 3 decimal places only.
# TODO: Investigate this.
