# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

#-------------------------------------------------- -
# all endpoints normal
getMVNorm <- function(nSubject = 100,
                      vNormMean = c(0.1, 0.4),
                      mNormSigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
                      nSeed = 123) {
  set.seed(seed = nSeed)
  mvtnorm::rmvnorm(
    n = nSubject,
    mean = as.numeric(vNormMean),
    sigma = as.matrix(mNormSigma)
  )
}


#------------------------------------------------- -
# all endpoints binary
getMVBinom <- function(nSubject = 100,
                       vProp = c(0.1, 0.2),
                       mNormCorr = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
                       nSeed = 123) {
  set.seed(seed = nSeed)
  mResNorm <- mvtnorm::rmvnorm(
    n = nSubject,
    sigma = as.matrix(mNormCorr)
  )

  vCutOffs <- qnorm(1 - vProp)
  t(apply(mResNorm, 1, function(x) {
    as.integer(x >= vCutOffs)
  }))
}

#------------------------------------------------- -
# endpoints either normal or binary
getMVNormBinom <- function(nSubject = 100000,
                           vProp = c(0.1, 0.2),
                           mNormCorr = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
                           nSeed = 123) {
  set.seed(seed = nSeed)
  mvtnorm::rmvnorm(
    n = nSubject,
    mean = as.numeric(vNormMean),
    sigma = as.matrix(mNormSigma)
  )
}
