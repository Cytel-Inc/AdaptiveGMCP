# FS_GMCP_Example.R
# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

#######################################################################
#####   Script to Run Examples for Fixed Sample GMCP Test  #############
#######################################################################

library(AdaptGMCP)

##############################################################
# Examples from the following paper:
# "A graphical approach to sequentially rejective multiple test procedures",
# Bretz et. al. 2009
# File: Bretz-Graphical-MCP- 2009b.pdf

# EXAMPLE 1 #################################
# From pages 4 and 5 of Bretz's paper:
# 3 hypotheses to test using the Bonferroni-Holm test and graphical method
# Equal initial allocation of significance level to the 3 hypo
# All hypotheses at the same level, i.e. no hierarchy among them
# Total type I error: 0.05

# Setting input parameters for the function

# Weights
wi <- c(1/3, 1/3, 1/3)

# Transition matrix
g <- matrix(c(0, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0), byrow = T, nrow = 3)

# Test type
test <- "Bonf" # Bonferroni

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
# Enter these raw p-values when asked: p1=0.02, p2=0.055, and p3=0.012
out <- adaptGMCP_PC(WI=wi, G=g, test.type=test, alpha = alp, info_frac = t,
                    plotGraphs = T)

# EXAMPLE 2 ###################################
# From page 11 of Bretz's paper:
# 3 hypotheses: H1, H2 are of primary interest and H3 is of interest only if
# H1 and H2 can be both rejected.
# Perform the Bonferroni–Holm procedure at level alpha for H1 and H2. If both
# H1 and H2 can be rejected, then H3 is tested at level alpha.
# Total alpha = 0.05
# Observed p-values are p1=0.04, p2=0.01, and p3=0.03.

# Setting input parameters for the function

# Weights
wi <- c(1/2, 1/2, 0)

# Transition matrix
eps <- 1e-6
g <- matrix(c(0, 1, 0, 1-eps, 0, eps, 0, 0, 0), byrow = T, nrow = 3)

# Test type
test <- "Bonf" # Bonferroni

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
# Enter these raw p-values when asked: p1=0.04, p2=0.01, and p3=0.03
out <- adaptGMCP_PC(WI=wi, G=g, test.type=test, alpha = alp, info_frac = t,
                    plotGraphs = T)

# EXAMPLE 3 ######################
# From page 12 of Bretz's paper
# 4 hypotheses: H3 and H4 are of interest only if both H1 and H2 were rejected.
# alpha=0.05, r1=0.8, r2=0.2
# Raw p-values: p1=0.04, p2=0.01, p3=0.03, and p4=0.04

# Setting input parameters for the function

# Weights
wi <- c(1/2, 1/2, 0, 0)

# Transition matrix
eps <- 1e-6
r1 <- 0.8
r2 <- 0.2
g <- matrix(c(0, 1, 0, 0,
              1-eps, 0, r1*eps, r2*eps,
              0, 0, 0, 1,
              0, 0, 1, 0), byrow = T, nrow = 4)

# Test type
test <- "Bonf" # Bonferroni

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
# Enter these raw p-values when asked: p1=0.04, p2=0.01, p3=0.03, and p4=0.04
out <- adaptGMCP_PC(WI=wi, G=g, test.type=test, alpha = alp, info_frac = t,
                    plotGraphs = T)

# EXAMPLE 4 ################################
# From page 13 of Bretz's paper
# First with the graph from Figure 2: gatekeeping procedure with 4 hypotheses
# H1 and H2 are primary hypotheses and H3 and H4 are secondary hypotheses.
# alpha = 0.05
# Raw p-values: p1=0.02, p2=0.04, p3=0.01, and p4=0.015

# Setting input parameters for the function

# Weights
wi <- c(0.5, 0.5, 0, 0)

# Transition matrix
g <- matrix(c(0, 0, 0.5, 0.5,
              0, 0, 0.5, 0.5,
              0, 0, 0, 1,
              0, 0, 1, 0), byrow = T, nrow = 4)

# Test type
test <- "Bonf" # Bonferroni

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
# Enter these raw p-values when asked: p1=0.02, p2=0.04, p3=0.01, and p4=0.015
out <- adaptGMCP_PC(WI=wi, G=g, test.type=test, alpha = alp, info_frac = t,
                    plotGraphs = T)

# Now with the improved gatekeeping procedure suggested by Bretz et al

# Transition matrix
eps <- 1e-6
g <- matrix(c(0, 0, 0.5, 0.5,
              0, 0, 0.5, 0.5,
              eps, 0, 0, 1-eps,
              0, eps, 1-eps, 0), byrow = T, nrow = 4)

# Calling the function
# Enter these raw p-values when asked: p1=0.02, p2=0.04, p3=0.01, and p4=0.015
out <- adaptGMCP_PC(WI=wi, G=g, test.type=test, alpha = alp, info_frac = t,
                    plotGraphs = T)

# EXAMPLE 5 ######################################
# Ref: Bretz's 2011 paper "Graphical approaches for multiple comparison
# procedures using weighted Bonferroni, Simes, or parametric tests"
# Pages 14-15
# Trials comparing a new compound with placebo for two primary and two secondary
# endpoints, i.e. 4 hypotheses
# H3 and H4 are tested only if at least one of H1 and H2 are rejected.
# H1 and H2 are tested with truncated Holm procedure and H3 and H4 with regular
# Holm procedure.
# Raw p-values: p1=0.0121. p2=0.0337, p3=0.0084, p4=0.0160
# Gamma=0.5 for Truncated Holm procedure.
# alpha=0.05

# Setting input parameters for the function

# Weights
wi <- c(0.5, 0.5, 0, 0)

# Transition matrix
gm <- 0.5 # Gamma for Truncated Holm
g <- matrix(c(0, gm, (1-gm)/2, (1-gm)/2,
              gm, 0, (1-gm)/2, (1-gm)/2,
              0, 0, 0, 1,
              0, 0, 1, 0), byrow = T, nrow = 4)

# Test type
test <- "Bonf" # Bonferroni

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
# Enter these raw p-values when asked: p1=0.0121, p2=0.0337, p3=0.0084, p4=0.0160
out <- adaptGMCP_PC(WI=wi, G=g, test.type=test, alpha = alp, info_frac = t,
                    plotGraphs = T)
######################################################

# EXAMPLE 6 ###########################################
# Ref: Bretz's 2011 paper "Graphical approaches for multiple comparison
# procedures using weighted Bonferroni, Simes, or parametric tests"
# Page 6 onwards
# 4 hypotheses
# Raw p-values: p1=0.01, p2=0.005, p3=0.1, p4=0.5
# alpha=0.025

# Setting input parameters for the function

# Weights
wi <- c(0.5, 0.5, 0, 0)

# Transition matrix
g <- matrix(c(0, 0, 1, 0,
              0, 0, 0, 1,
              0, 1, 0, 0,
              1, 0, 0, 0), byrow = T, nrow = 4)

# Test type
test <- "Bonf" # Bonferroni

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
# Enter these raw p-values when asked: p1=0.01, p2=0.005, p3=0.1, p4=0.5
out <- adaptGMCP_PC(WI=wi, G=g, test.type=test, alpha = alp, info_frac = t,
                    plotGraphs = T)

# EXAMPLE 7 #################################
# Ref: Bretz's 2011 paper "Graphical approaches for multiple comparison
# procedures using weighted Bonferroni, Simes, or parametric tests"
# Page 11, example 3
# 4 hypotheses: H1, H2 non-inf for low and high dose. H3, H4 sup for the same
# alpha = 0.025
# Correlations: H1xH3=1, H2xH4=1, all other = 0.5

# Setting input parameters for the function

# Weights
wi <- c(0.5, 0.5, 0, 0)

# Transition matrix
g <- matrix(c(0, 0, 1, 0,
              0, 0, 0, 1,
              0, 1, 0, 0,
              1, 0, 0, 0), byrow = T, nrow = 4)

# Test type
test <- "Dunnett"

# Type I error
alp <- 0.025

# Info fraction
t <- c(1) # This is a fixed sample design

corr <- matrix(c(1, 0.5, 1, 0.5,
                 0.5, 1, 0.5, 1,
                 1, 0.5, 1, 0.5,
                 0.5, 1, 0.5, 1), byrow = T, nrow = 4)

# We will ignore typeOfDesign argument since that is required in case of
# group sequential designs.
# We will also ignore other related arguments like deltaWT, deltaPT1, gammaA,
# and userAlphaSpending.

# Similarly, we will ignore other parameters like MultipleWinners,
# Selection, and UpdateStrategy.

# We will set plotGraphs to TRUE to produce the graphs.

# Calling the function
# Enter these raw p-values when asked: p1=0.01, p2=0.02, p3=0.005, p4=0.5
out <- adaptGMCP_PC(WI=wi, G=g, test.type=test, alpha = alp, info_frac = t,
                    Correlation = corr, plotGraphs = T)

# EXAMPLE 8 ######################################
# Ref: Bretz's 2011 paper "Graphical approaches for multiple comparison
# procedures using weighted Bonferroni, Simes, or parametric tests"
# Page 14, last paragraph from the Weighted Simes Test section
# 4 hypotheses
# Raw p-values: p1=0.01, p2=0.005, p3=0.015, p4=0.022
# alpha=0.025

# Setting input parameters for the function

# Weights
wi <- c(0.5, 0.5, 0, 0)

# Transition matrix
g <- matrix(c(0, 0, 1, 0,
              0, 0, 0, 1,
              0, 1, 0, 0,
              1, 0, 0, 0), byrow = T, nrow = 4)

# Test type
test <- "Bonf" # Bonferroni

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
# Enter these raw p-values when asked: p1=0.01, p2=0.005, p3=0.015, p4=0.022
out <- adaptGMCP_PC(WI=wi, G=g, test.type=test, alpha = alp, info_frac = t,
                    plotGraphs = T)

# Now trying the same example with Simes test

# Test type
test <- "Simes"

# Calling the function
# Enter these raw p-values when asked: p1=0.01, p2=0.005, p3=0.015, p4=0.022
out <- adaptGMCP_PC(WI=wi, G=g, test.type=test, alpha = alp, info_frac = t,
                    plotGraphs = T)

