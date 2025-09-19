# try_gmcp.R
# Experimenting with the package gMCP in this file

# library(dplyr)
library(gMCP)
# library(bindata)
# library(purrr)
# library(mvtnorm)

# # BINOMIAL SIMULATIONS ######################################
# Power simulations for continuous response
mGraphMat <- matrix(rep(0, 16), byrow = T, nrow = 4)
colnames(mGraphMat) <- rownames(mGraphMat) <- c("H1", "H2", "H3", "H4")

gGraph <- new("graphMCP", m=mGraphMat, weights=rep(1/4, 4))
gGraph

dArmProps <- c(0.1, 0.1, 0.1, 0.1, 0.1)
nArmSS <- rep(100, 5)
z_means <- (dArmProps[-1] - dArmProps[1]) /
  sqrt(dArmProps[-1]*(1-dArmProps[-1])/nArmSS[-1] +
         dArmProps[1]*(1-dArmProps[1])/nArmSS[1])
z_means

mCorr <- matrix(0.5, nrow = 4, ncol = 4)
diag(mCorr) <- 1
mCorr

dAlpha <- 0.025

nSim <- 100 # 10000

# Running the simulations
# Using parametric graphical MCP
system.time({
  out1 <- calcPower(graph=gGraph, alpha=dAlpha, mean=z_means, corr.sim = mCorr,
                    corr.test = mCorr, n.sim = nSim)
})

#10k simulations with parametric option took 412 seconds (about 6.5 mins).
out1

# Running the simulations
# Using parametric graphical MCP
system.time({
  out2 <- calcPower(graph=gGraph, alpha=dAlpha, mean=z_means, corr.sim = mCorr,
                    corr.test = NULL, n.sim = nSim)
})

out2

# # NOTE: Use the following method if the underlying distribution of the test
# # statistics is NOT multivariate normal.
# # E.g. the following example is for binary endpoints with relatively small
# # sample size of 20 per arm. The resulting binomial distribution is not well
# # approximated by the normal distribution.
#
# # We have 2 endpoints and 2 treatment arms (and one common control arm).
# # Thus, we have 4 hypotheses:
# # H1: EP1, Trtm1 vs Ctrl
# # H2: EP1, Trtm2 vs Ctrl
# # H3: EP2, Trtm1 vs Ctrl
# # H4: EP2, Trtm2 vs Ctrl
#
# # Correlation matrix for the resulting multivariate binomial distribution
# # corresponding to these 4 hypo:
# mCr <- matrix(c(1, 0.5, 0.3, 0.15,
#                0.5, 1, 0.15, 0.3,
#                0.3, 0.15, 1, 0.5,
#                0.15, 0.3, 0.5, 1), byrow = T, nrow = 4)
#
# # Alternative distributions to be used for simulating data
# dProps <- c(0.35, 0.4, 0.25, 0.3)
#
# n1 <- 20 # binomial n
#
# # # Generating a sample from the multivariate binomial distribution
# # # Return value is a matrix with one column per hypothesis, i.e. it is a n1x4 matrix.
# # mBinVec <- rmvbin(n1, margprob = dProps, bincorr = mCr)
# # mBinVec
# #
# # # Binomial sample sums
# # dSmplSums <- colSums(mBinVec)
# # dSmplSums
# #
# # # Calculating raw/marginal p-values for the 4 hypo using binom.test()
# # lTestOut <- map(dSmplSums, binom.test, n1, alternative="less")
# # lTestOut
# #
# # dRawPVal <- map_dbl(lTestOut, function(x) x$p.value)
# # dRawPVal
#
# # Now writing a function for doing the above steps
# CalcRawPVals <- function(nBin_n, dBin_Props, mBinCorr, alt){
#   # Generating a sample from the multivariate binomial distribution
#   # Return value is a matrix with one column per hypothesis, i.e. it is a n1x4 matrix.
#   mBinVec <- rmvbin(nBin_n, margprob = dBin_Props, bincorr = mBinCorr)
#
#   # Binomial sample sums
#   dSmplSums <- colSums(mBinVec)
#
#   # Calculating raw/marginal p-values for the 4 hypo using binom.test()
#   lTestOut <- map(dSmplSums, binom.test, nBin_n, alternative=alt) # "less")
#
#   dRawPVal <- map_dbl(lTestOut, function(x) x$p.value)
#
#   return(dRawPVal)
# }
#
# nSims <- 10000
# set.seed(123)
# system.time({
#   dPVals <- replicate(nSims, CalcRawPVals(n1, dProps, mCr, "less"))
# })
# dPVals
# dPVals <- t(dPVals)
# dPVals
#
# # nSims <- 1000 # 10000
# # system.time({
# #   pvals <- t(
# #   replicate(nSims,
# #             sapply(colSums(rmvbin(n1, margprob = dProps, bincorr = mCr)),
# #                    function(x, ...) {binom.test(x, ...)$p.value},
# #                    n=n1, alternative="less")
# # ))})
# #
# # pvals
#
# lGraph <- generalSuccessive(gamma = 0, delta = 0)
# lGraph
#
# lOut <- graphTest(pvalues = dPVals, graph = lGraph)
# lOut
# extractPower(lOut)
#
# # BINOMIAL SIMULATIONS OVER >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#

# # NORMAL SIMULATIONS ##############################
# # Here it is assumed that the distribution of the test statistics is multivariate
# # normal. We also make this assumption in AdaptGMCP.
#
# ## reproduce example from Stat Med paper (Bretz et al. 2010, Table I)
# ## first only consider line 2 of Table I
# ## significance levels
# graph <- simpleSuccessiveII()
# ## alternative (mvn distribution)
# corMat <- rbind(c(1, 0.5, 0.5, 0.5/2),
#                 c(0.5,1,0.5/2,0.5),
#                 c(0.5,0.5/2,1,0.5),
#                 c(0.5/2,0.5,0.5,1))
# theta <- c(3, 0, 0, 0)
# system.time({
#   out <- calcPower(graph=graph, alpha=0.025, mean=theta, corr.sim=corMat, n.sim= 10000) # 100000)
# })
# out
#
# # Try_mvtnorm <- function(nArms, nPerArmSS, nSims, mu, sig){
# #   for(i in 1:nSims){
# #     for(j in 1:nArms){
# #       rmvnorm(nPerArmSS, mean=mu, sigma = sig)
# #     }
# #   }
# # }
# #
# # system.time({
# #   Try_mvtnorm(nArms = 4, nPerArmSS = 100, nSims = 100000, mu=theta, sig = corMat)
# # })

# LEARNING TO USE GMCP ######################################
# Create a new graph
graph <- matrix(c(0, 0.5, 0.5,
                  0, 0, 0,
                  0, 0, 0), byrow = T, nrow = 3)
rownames(graph) <- colnames(graph) <- c("H1", "H2", "H3")

# Now creating the graph
g <- new("graphMCP", m=graph, weights=c(1, 0, 0))
g

# Run gMCP procedure with data
z <- c(2.8, 1.5, 2.1)
# Correlation matrix
corMat <- matrix(rep(0.5, 9), nrow = 3, byrow = T)
diag(corMat) <- 1

result1 <- gMCP(g, pvalues = 2*(1-pnorm(abs(z))), correlation = corMat, verbose = T)
result1@adjPValues

result2 <- gMCP(g, pvalues = 2*(1-pnorm(abs(z))), verbose = T)
result2@adjPValues

# Trying with correlated test statistics
# NOTE: when a correlation matrix is specified to gMCP, it is assumed that the
# Z test stats follow a multivariate normal distribution with the given
# correlation matrix.
# If p1, p2, ... are the raw p-values, then the corresponding Z test stats are
# phiinv(1-p1), phiinv(1-p2), ...
# Note that gMCP supports partial knowledge of the correlation matrix. So, some
# matrix entries can be specified and some can be set to NA. E.g. this is the
# case when you are comparing 2 doses of a treatment to common control and are
# observing one primary and one second endpoint. Here H1 and H2 are hypo for the
# high dose and low dose (both for the primary endpoint), whereas H3 and H4 are
# for the high dose and low dose (both for the secondary endpoint). In this case,
# the correlation between H1 and H2 and that between H3 and H4 are known, but
# that between the other hypothesis pairs is unknown.

corr <- matrix(c(1, 0.5, 0.3,
                 0.5, 1, 0.4,
                 0.3, 0.4, 1), byrow = T, nrow = 3) # Correlation matrix for test stats

result3 <- gMCP(g, pvalues = 2*(1-pnorm(abs(z))), correlation = corr)
result3

result1@adjPValues
result2@adjPValues
result3@adjPValues

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Example: Testing multiple trtms vs a common control, with continuous outcomes
graphMat <- matrix(c(0, 0.5, 0.5,
                     0.5, 0, 0.5,
                     0.5, 0.5, 0), byrow = 3, nrow = 3)
colnames(graphMat) <- rownames(graphMat) <- c("H1", "H2", "H3")

# Creating the graph
gr <- new("graphMCP", m=graphMat, weights=rep(1/3,3))
gr

# z stats
z <- c(2.6, 2.2, 1.5)
p <- 2*(1-pnorm(abs(z)))
p

# Correlation matrix for the test stats
# The test stats are correlated due to the shared common control
# Under balanced randomization, the correlation between each hypotheses pair
# is 0.5.
corrMat <- matrix(c(1, 0.5, 0.5,
                    0.5, 1, 0.5,
                    0.5, 0.5, 1), byrow = T, nrow = 3)
colnames(corrMat) <- rownames(corrMat) <- c("H1", "H2", "H3")
corrMat

result4 <- gMCP(graph = gr, pvalues = p, correlation = corrMat)
result4

result5 <- gMCP(graph = gr, pvalues = p, test = "parametric", correlation = corrMat)
result5

result6 <- gMCP(graph = gr, pvalues = p, test = "Simes")
result6

result7 <- gMCP(graph = gr, pvalues = p, test = "Bonferroni")
result7

# Binomial example >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
x <- c(40, 60, 55, 50)
n <- rep(100, 4)
p <- x/n

# We want to test H1, H2, H3 where Hi:pi=p0.

# Parametric test
# First we need to calculate the z stats and then p-values.
z <- (p[-1] - p[1]) / sqrt(p[-1]*(1-p[-1])/n[-1] + p[1]*(1-p[1])/n[1])
pval <- 2*(1-pnorm(abs(z)))
#
# pooled <- (x[-1]+x[1])/(n[-1]+n[1])
# pooled <- c(0, pooled)
#
# z <- (p[-1] - p[1]) / sqrt(pooled[-1]*(1-pooled[-1])*(1/n[-1] + 1/n[1]))
# pval <- 2*(1-pnorm(abs(z)))

names(z) <- names(pval) <- c("H1", "H2", "H3")

# Correlation matrix for the test stats
# The test stats are correlated due to the shared common control
# Under balanced randomization, the correlation between each hypotheses pair
# is 0.5.
corrMat <- matrix(c(1, 0.5, 0.5,
                    0.5, 1, 0.5,
                    0.5, 0.5, 1), byrow = T, nrow = 3)
colnames(corrMat) <- rownames(corrMat) <- c("H1", "H2", "H3")
corrMat

graphMat <- matrix(c(0, 0.5, 0.5,
                     0.5, 0, 0.5,
                     0.5, 0.5, 0), byrow = 3, nrow = 3)
colnames(graphMat) <- rownames(graphMat) <- c("H1", "H2", "H3")

# Creating the graph
gr <- new("graphMCP", m=graphMat, weights=rep(1/3,3))
gr

# Performing a parametric test
result8 <- gMCP(graph = gr, pvalues = pval, test = "parametric",
                correlation = corrMat)
result8
result8@adjPValues

# Performing a non-parametric test
result9 <- gMCP(graph = gr, pvalues = pval, test = "Bonferroni")
result9
result9@adjPValues

# Performing Simes test
result10 <- gMCP(graph = gr, pvalues = pval, test = "Simes")
result10
result10@adjPValues

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Power simulations for continuous response
graphMat <- matrix(c(0, 0.5, 0.5,
                     0.5, 0, 0.5,
                     0.5, 0.5, 0), byrow = 3, nrow = 3)
colnames(graphMat) <- rownames(graphMat) <- c("H1", "H2", "H3")

graph <- new("graphMCP", m=graphMat, weights=rep(1/3, 3))

means <- c(3.0, 2.0, 0.0)

corr <- matrix(0.5, nrow = 3, ncol = 3)
diag(corr) <- 1

# Running the simulations
# Using parametric graphical MCP
system.time({
  out1 <- calcPower(graph=graph, alpha=0.05, mean=means, corr.sim = corr,
                    corr.test = corr, n.sim = 10000)
})

#10k simulations with parametric option took 412 seconds (about 6.5 mins).
out1

# Using non-parametric graphical MCP
# If argument corr.test is not specified, gMCP uses Bonferroni.
system.time({
  out2 <- calcPower(graph=graph, alpha=0.05, mean=means, corr.sim = corr,
                    n.sim = 10000)
})

# 10k simulations with non-parametric option took just 2 seconds.
out2

out1$LocalPower
out2$LocalPower
