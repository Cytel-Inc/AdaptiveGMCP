# Compute adjusted p-value for a parametric test
# p2 stage-2 p-values for the intersection hypothesis
# stage-1 boundary for the intersection hypothesis
# w2 stage-2 weights for the intersection hypothesis
# SigmaZ variance covariance matrix on Z scale for the intersection hypothesis
getParamAdjPValue <- function(
    p2 = c(0.02, 0.1),
    a1 = c(0.000782, 0.000782),
    w2 = c(0.5, 0.5),
    SigmaZ = matrix(
      c(
        1, 0.5, 0.7018299, 0.3535534,
        0.5, 1, 0.3535534, 0.7018299,
        0.7018299, 0.3535534, 1, 0.5,
        0.3535534, 0.7018299, 0.5, 1
      ),
      nrow = 4, byrow = T
    )) {
  q <- min(p2 / w2)
  upper <- c(qnorm(1 - a1), qnorm(1 - w2 * q))
  lower <- -Inf
  prob <- mvtnorm::pmvnorm(
    lower = lower,
    upper = upper,
    sigma = SigmaZ
  )[1]
  (1 - prob) # Under null this should be cummulative alpha for that look
}
