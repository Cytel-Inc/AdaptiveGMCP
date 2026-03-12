# Example for analysis in the population enrichment case
library(AdaptGMCP)

# EXAMPLE 1 #############################################
# Problem: Full population and a subpopulation (50% of full population)
# High dose and low dose of the drug being tested
# 2 stage trial with interim look at 50% information fraction
# H1: full pop, high dose; H2: subpop, high dose;
# H3: full pop, low dose; H4: subpop, low dose
# Correlations:
#       (1) For the same dose, test stats for full population and subpopulation
#           are correlated due to the shared patients and the correlation coefficient
#           is the same as the subpop prop. (H1 vs H2 and H3 vs H4)
#       (2) For the same population, test stats for high dose and low dose are correlated
#           due to the shared control arm and the correlation coefficient depends on
#           the allocation ratio between the control and the dose arms.
#           (H1 vs H3 and H2 vs H4)
#       (3) Correlation unknown between the remaining test stats (i.e. H1 vs H4 and H2 vs H3)

# Setting input parameters for the function

# Weights
wi <- rep(0.25, 4) # Initial weights for the 4 hypo

# Transition matrix
g <- matrix(c(0, 0.5, 0.5, 0, # H1->H2, H1->H3
              0.5, 0, 0, 0.5, # H2->H1, H2->H4
              0.5, 0, 0, 0.5, # H3->H1, H3->H4
              0, 0.5, 0.5, 0), # H4->H2, H4->H3
              byrow = T, nrow = 4)

# Test type
test <- "Dunnett" # "Bonf" # "Partly-Parametric" #

# Type I error
alp <- 0.025

# Info fraction
t <- c(0.5, 1) # This is a fixed sample design

# Design type
des <- "asOF"

# Correlation matrix between test stats
# corr <- matrix(c(1, 0.5, 0.5, NA,
#                  0.5, 1, NA, 0.5,
#                  0.5, NA, 1, 0.5,
#                  NA, 0.5, 0.5, 1), byrow = T, nrow = 4)
corr <- matrix(c(1, NA, 0.5, NA,
                 NA, 1, NA, 0.5,
                 0.5, NA, 1, NA,
                 NA, 0.5, NA, 1), byrow = T, nrow = 4)

# Calling the function
# Use these p-values: p(H1) = 0.00025, p(H2) = 0.0952, p(H3) = 0.0245, p(H4) = 0.1104
out <- adaptGMCP_PC(WI=wi, G=g, test.type = test, alpha = alp, info_frac = t,
                    typeOfDesign = des, Correlation = corr,
                    Selection = T, UpdateStrategy = T, plotGraphs = T)

print(out)
