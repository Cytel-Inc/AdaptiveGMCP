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
corr <- matrix(c(1, 0.5, 0.5, NA,
                 0.5, 1, NA, 0.5,
                 0.5, NA, 1, 0.5,
                 NA, 0.5, 0.5, 1), byrow = T, nrow = 4)
# corr <- matrix(c(1, NA, 0.5, NA,
#                  NA, 1, NA, 0.5,
#                  0.5, NA, 1, NA,
#                  NA, 0.5, NA, 1), byrow = T, nrow = 4)
x <- conn.comp(corr)
print(x)

y <- clique.partition(corr)
print(y)

# Calling the function
# Use these p-values: p(H1) = 0.00025, p(H2) = 0.0952, p(H3) = 0.0245, p(H4) = 0.1104
out <- adaptGMCP_PC(WI=wi, G=g, test.type = test, alpha = alp, info_frac = t,
                    typeOfDesign = des, Correlation = corr,
                    Selection = T, UpdateStrategy = T, plotGraphs = T)

# EXAMPLE 2 #############################################
# 8-hypothesis population enrichment problem: 2 doses (dose 1, dose 2) x 2 endpoints
# (primary, secondary) x 2 populations (full population, HPV+ subgroup).
# H1/H2: primary endpoint, full population, dose 1/2 (initial weight 0.35 each).
# H3/H4: secondary endpoint, full population, dose 1/2 (initial weight 0).
# H5/H6: primary endpoint, HPV+ subgroup, dose 1/2 (initial weight 0.15 each).
# H7/H8: secondary endpoint, HPV+ subgroup, dose 1/2 (initial weight 0).
G <- matrix(
  c(
  # H1    H2    H3    H4    H5    H6    H7    H8
    0,    0.2,  0.4,  0,    0.2,  0.2,  0,    0,    # H1
    0.2,  0,    0,    0.4,  0.2,  0.2,  0,    0,    # H2
    0,    1/3,  0,    0,    1/3,  1/3,  0,    0,    # H3
    1/3,  0,    0,    0,    1/3,  1/3,  0,    0,    # H4
    0.2,  0.2,  0,    0,    0,    0.2,  0.4,  0,    # H5
    0.2,  0.2,  0,    0,    0.2,  0,    0,    0.4,  # H6
    1/3,  1/3,  0,    0,    0,    1/3,  0,    0,    # H7
    1/3,  1/3,  0,    0,    1/3,  0,    0,    0     # H8
  ),
  nrow = 8, byrow = TRUE,
  dimnames = list(
    c("H1","H2","H3","H4","H5","H6","H7","H8"),
    c("H1","H2","H3","H4","H5","H6","H7","H8")
  )
)

corr <- matrix(
  c(
  # H1    H2    H3    H4    H5    H6    H7    H8
    1,    0.5,  NA,   NA,   1,    0.5,  NA,   NA,   # H1
    0.5,  1,    NA,   NA,   0.5,  1,    NA,   NA,   # H2
    NA,   NA,   1,    0.5,  NA,   NA,   1,    0.5,  # H3
    NA,   NA,   0.5,  1,    NA,   NA,   0.5,  1,    # H4
    1,    0.5,  NA,   NA,   1,    0.5,  NA,   NA,   # H5
    0.5,  1,    NA,   NA,   0.5,  1,    NA,   NA,   # H6
    NA,   NA,   1,    0.5,  NA,   NA,   1,    0.5,  # H7
    NA,   NA,   0.5,  1,    NA,   NA,   0.5,  1     # H8
  ),
  nrow = 8, byrow = TRUE,
  dimnames = list(
    paste0("H", 1:8),
    paste0("H", 1:8)
  )
)
