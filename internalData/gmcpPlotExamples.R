# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

# ## Example-1: Primary and secondary type hypothesis
# #Initial weights
# WI <- c(1/2,1/2,0,0)
#
# #Transition Matrix
# G <- matrix(c(0,0.75,0.25,0,
#               0.75,0,0,0.25,
#               0,1,0,0,
#               1,0,0,0),
#             nrow = 4, byrow = TRUE)
# #Group of hypothesis
# hGroup <- c('P','P','S','S')
#
# #Coordinates of the nodes(Hypothesis)
# cordinates <- list(c(-1,1),c(1,1),c(-1,-1),c(1,-1))
#
# #Plot
# gmcpPlot(WI = WI, G = G, hGroup = hGroup, cordinates = cordinates)
#
#
# ##################################################################
# ## Example-2: Sequential rejection strategy for multi-arm study
#
# #Initial weights
# WI <- c(1/5,1/5,1/5,1/5,1/5)
#
# #Transition Matrix
# G  <- matrix(c(0,1/4,1/4,1/4,1/4,
#                1/4,0,1/4,1/4,1/4,
#                1/4,1/4,0,1/4,1/4,
#                1/4,1/4,1/4,0,1/4,
#                1/4,1/4,1/4,1/4,0),
#              nrow = 5, byrow = TRUE)
#
# #Plot
# gmcpPlot(WI = WI, G = G)
#
