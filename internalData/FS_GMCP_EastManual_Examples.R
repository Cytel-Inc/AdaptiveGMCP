# File: FS_GMCP_EastManual_Examples.R
# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

#######################################################################
#####   Script to Run Examples for Fixed Sample GMCP Test  #############
#######################################################################

library(AdaptGMCP)

# We will now try examples from the East user manual.
# Ref: Chapter 100: Analysis‐Multiple Comparison Procedures for Continuous Data

# EXAMPLE M1 ##################################
# From page 2550 of the manual
# 4 treatment arms vs control
# No hierarchy among the 4 hypotheses (this chapter is about plain fixed sample
# MCP tests)
# Test type: Dunnett
# alpha = 0.025
# Raw p-values: p1=0.638, p2=0.01, p3=0.007, p4=3.959E-4

# Setting input parameters for the function

# Weights
wi <- c(1/4, 1/4, 1/4, 1/4)

# Transition matrix
g <- matrix(c(0, 1/3, 1/3, 1/3,
              1/3, 0, 1/3, 1/3,
              1/3, 1/3, 0, 1/3,
              1/3, 1/3, 1/3, 0), byrow = T, nrow = 4)

# Test type
test <- "Dunnett"

# Type I error
alp <- 0.025

# Info fraction
t <- c(1) # This is a fixed sample design

corr <- matrix(c(1, 0.5, 0.5, 0.5,
                 0.5, 1, 0.5, 0.5,
                 0.5, 0.5, 1, 0.5,
                 0.5, 0.5, 0.5, 1), byrow = T, nrow = 4)

# We will ignore typeOfDesign argument since that is required in case of
# group sequential designs.
# We will also ignore other related arguments like deltaWT, deltaPT1, gammaA,
# and userAlphaSpending.

# Similarly, we will ignore other parameters like Correlation, MultipleWinners,
# Selection, and UpdateStrategy.

# We will set plotGraphs to TRUE to produce the graphs.

# Calling the function
# Enter these raw p-values when asked: p1=0.638, p2=0.01, p3=0.007, p4=3.959E-4
out <- adaptGMCP_PC(WI=wi, G=g, test.type=test, alpha = alp, info_frac = t,
                    Correlation = corr, plotGraphs = T)
# out <- adaptGMCP_PC(WI=wi, G=g, test.type=test, alpha = alp, info_frac = t,
#                     plotGraphs = T)

# EXAMPLE M2 #######################
# From chapter 101, pages 2579 onwards
# 4 doses compared against placebo, no hierarchy among them
# 2 scenarios: increasing DR profile, decreasing DR profile
# Scenario 1 raw p-values from table 101.57:
#           p1=0.638138, p2=0.009838, p3=0.00673, p4=0.000396
# alpha=0.025

# Setting input parameters for the function

# Weights
wi <- c(1/4, 1/4, 1/4, 1/4)

# Transition matrix
g <- matrix(c(0, 1/3, 1/3, 1/3,
              1/3, 0, 1/3, 1/3,
              1/3, 1/3, 0, 1/3,
              1/3, 1/3, 1/3, 0), byrow = T, nrow = 4)

# We will try different test types and compare output with table 101.59.
# Test type
test <- "Dunnett"

# Type I error
alp <- 0.025

# Info fraction
t <- c(1) # This is a fixed sample design

# We will ignore typeOfDesign argument since that is required in case of
# group sequential designs.
# We will also ignore other related arguments like deltaWT, deltaPT1, gammaA,
# and userAlphaSpending.

# Similarly, we will ignore other parameters like Correlation, MultipleWinners,
# Selection, and UpdateStrategy.

# We will set plotGraphs to TRUE to produce the graphs.

# Calling the function
# Enter these raw p-values when asked:
#         p1=0.638138, p2=0.009838, p3=0.00673, p4=0.000396
out <- adaptGMCP_PC(WI=wi, G=g, test.type=test, alpha = alp, info_frac = t,
                    plotGraphs = T)

# WARNING: When we compare the output of above with that from the table 101.59,
# We find that only the adjusted p-value for H1 matches and that too under
# Step down Dunnett test in the table.

# Trying out with Bonferroni test:
test <- "Bonf"

# Calling the function
# Enter these raw p-values when asked:
#         p1=0.638138, p2=0.009838, p3=0.00673, p4=0.000396
out <- adaptGMCP_PC(WI=wi, G=g, test.type=test, alpha = alp, info_frac = t,
                    plotGraphs = T)

# WARNING: The output of the above is not matching with what the user manual
# mentions. I need to investigate that and then try the other tests mentioned
# in that section of the manual.

#############################
# Checking how the gMCP package handles this problem...
# TODO: TRY THIS AFTER INSTALLING JAVA.
#
# # Install gMCP if not already installed
# # install.packages("gMCP")
#
# library(gMCP)
#
# # Define the transition matrix: 4x4 with 0.33 between hypotheses (0 on diagonals)
# transition.matrix <- matrix(
#   c(0,   1/3, 1/3, 1/3,
#     1/3, 0,   1/3, 1/3,
#     1/3, 1/3, 0,   1/3,
#     1/3, 1/3, 1/3, 0),
#   nrow = 4, byrow = TRUE
# )
#
# # Define weights for each hypothesis (equal: 0.25 each)
# weights <- rep(0.25, 4)
#
# # Create the graph
# graph <- new("GraphMCP", m=transition.matrix, weights=weights)
#
# # Plot the graph
# plot(graph, main = "Graphical Bonferroni Procedure (Equal Weights)")

#############################

# EXAMPLE M3 #####################################
# Now trying tests for multiple endpoints from the East user manual.
# Reference: East user manual > chapter 102 > pages 2586 onwards (starting with
# the screenshot on page 2586).
# Comparison of one treatment against placebo.
# However, there are 2 primary (H1, H2) and 2 secondary (H3, H4) hypotheses.
# In this case, the hypotheses are typically hierarchical.
# Serial gatekeeping: If both H1 and H2 are rejected, then only H3 and H4 tested.
# Type I error: 0.05
# Raw p-values: p1=0.076, p2=0.035, p3=0.563, p4=0.407

# Setting input parameters for the function

# Weights
wi <- c(1/2, 1/2, 0, 0)

# Transition matrix
eps <- 1e-6
g <- matrix(c(0, 1-eps, eps/2, eps/2,
              1-eps, 0, eps/2, eps/2,
              0, 0, 0, 1,
              0, 0, 1, 0), byrow = T, nrow = 4)

# We will try different test types and compare output with table 101.59.
# Test type
test <- "Bonf"

# Type I error
alp <- 0.05

# Info fraction
t <- c(1) # This is a fixed sample design

# We will ignore typeOfDesign argument since that is required in case of
# group sequential designs.
# We will also ignore other related arguments like deltaWT, deltaPT1, gammaA,
# and userAlphaSpending.

# Similarly, we will ignore other parameters like Correlation, MultipleWinners,
# Selection, and UpdateStrategy.

# We will set plotGraphs to TRUE to produce the graphs.

# Calling the function
# Enter these raw p-values when asked:
#         p1=0.076, p2=0.035, p3=0.563, p4=0.407
out <- adaptGMCP_PC(WI=wi, G=g, test.type=test, alpha = alp, info_frac = t,
                    plotGraphs = T)

# WARNING: Results do not match with those from the manual. AdaptGMCP does not
# appear to work for serial gatekeeping.

