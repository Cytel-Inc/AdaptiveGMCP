# This function generate single arm(having multiple endpoint) responses given the normal scale correlation matrix
# nSubject number of response to generate
# lEpType list with endpoint types
# mNormCorr correlation matrix in normal scale(when marginals are not normal)
# mNormSigma covariance matrix in normal scale(when marginals are all normal)
# vNormMean vector of normal means(when marginals are all normal)
# vProp vector of proportions(when marginals are all binary)
# nSeed seed
# return: array of responses(columns = endpoints, Rows = subject)
genNormToOther <- function(nSubject = 100,
                           lEpType = list("Continuous", "Continuous"),
                           mNormCorr = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
                           mNormSigma = matrix(c(1, 0.5, 0.5, 1), nrow = 2),
                           vNormMean = c(0.1, 0.4),
                           vProp = c(0.1, 0.2),
                           nSeed = 123) {
  if (all(lEpType == "Continuous")) { # all endpoints continuous
    mResponse <- getMVNorm(
      nSubject = nSubject,
      vNormMean = vNormMean,
      mNormSigma = mNormSigma,
      nSeed = nSeed
    )
  } else if (all(lEpType == "Binary")) { # all endpoints binary
    mResponse <- getMVBinom(
      nSubject = nSubject,
      vProp = vProp,
      mNormCorr = mNormCorr,
      nSeed = nSeed
    )
  }
  # do Something
  mResponse
}


#------------------------------------------------- -
# This function generate single arm(having multiple endpoint) responses given the normal scale correlation matrix
# nArmID Arm Index
# nSubject number of response to generate
# vEPs Endpoints associated with the Arm Index(Existing)
# vEPType Endpoints type(Note Same length as vEPs)
# lNormMean User input continuous arm means
# lNormStdDev User input continuous arm Standard Deviation
# lProp User input binary arm proportions
# mNormCorr User input correlation matrix(normal scale)
# nSeed seed
genNormToOther2 <- function(nArmID = 3,
                            nSubject = 100,
                            vEPs = c(1, 4),
                            vEPType = c("Continuous", "Binary"),
                            lNormMean = list(
                              "EP1" = c(0, 0.4, 0.3),
                              "EP2" = NA,
                              "Ep3" = c(0.2, 0.55, 0.35),
                              "EP4" = NA
                            ),
                            lNormStdDev = list(
                              "EP1" = c(1.1, 1.2, 1.3),
                              "EP2" = NA,
                              "Ep3" = c(1.15, 1.25, 1.35),
                              "EP4" = NA
                            ),
                            lProp = list(
                              "EP1" = NA,
                              "EP2" = c(0.1, 0.2, 0.3),
                              "EP3" = NA,
                              "EP4" = c(0.5, 0.6, 0.7)
                            ),
                            mNormCorr = matrix(
                              c(
                                1, 0.5, 0.5, 0.5,
                                0.5, 1, 0.5, 0.5,
                                0.5, 0.5, 1, 0.5,
                                0.5, 0.5, 0.5, 1
                              ),
                              nrow = 4
                            ),
                            nSeed = 123) {
  # Mean for the ArmID(Continuous = user given, Binary = 0)
  vArmMean <- sapply(1:length(vEPs), function(epIDX) {
    if (vEPType[epIDX] == "Continuous") {
      lNormMean[[vEPs[epIDX]]][nArmID]
    } else if (vEPType[epIDX] == "Binary") {
      0
    }
  })

  #--------------------------------------------------------
  if(length(vEPs) == 1){
    #Single endpoint case
    if(vEPType[vEPs] == "Continuous"){
      mArmSigma <- lNormStdDev[[vEPs]][nArmID]^2
    }else if(vEPType[vEPs] == "Binary"){
      mArmSigma <- 1
    }
  }else{
    #Multiple endpoint case
    # Covariance Matrix for the ArmID
    # Continuous: covarince = sx*sy*rxy, Varince = sx*sx
    # Binary : covariance = rxy, Variance = 1
    mArmSigma <- matrix(NA, nrow = length(vEPs), ncol = length(vEPs))
    for (rowIDX in 1:length(vEPs)) {
      for (colIDX in rowIDX:length(vEPs)) {
        if (all(c(vEPType[rowIDX], vEPType[colIDX]) == "Continuous")) {
          #Continuous
          mArmSigma[colIDX, rowIDX] <- mArmSigma[rowIDX, colIDX] <-
            mNormCorr[vEPs[rowIDX], vEPs[colIDX]] *
            lNormStdDev[[vEPs[rowIDX]]][nArmID] *
            lNormStdDev[[vEPs[colIDX]]][nArmID]
        } else if (all(c(vEPType[rowIDX], vEPType[colIDX]) == "Binary")) {
          mArmSigma[colIDX, rowIDX] <- mArmSigma[rowIDX, colIDX] <-
            mNormCorr[vEPs[rowIDX], vEPs[colIDX]]
        }
      }
    }
  }
  set.seed(seed = nSeed)

  mNormRes <- mvtnorm::rmvnorm(
    n = nSubject,
    mean = as.numeric(vArmMean),
    sigma = as.matrix(mArmSigma)
  )

  if (all(vEPType == "Continuous")) {
    out <- mNormRes
  } else if (any(vEPType == "Binary")) {
    out <- NormToBin(
      nArmID = nArmID,
      vEPs = vEPs,
      vEPType = vEPType,
      lProp = lProp,
      vCutOffs = vCutOffs,
      mNormRes = mNormRes
    )
  }
  out
}

NormToBin <- function(nArmID, vEPs, vEPType, lProp, vCutOffs, mNormRes) {
  mMixRes <- mNormRes
  binIDX <- which(vEPType == "Binary")
  vProp <- sapply(vEPs[binIDX], function(epIDX) {
    lProp[[epIDX]][nArmID]
  })

  vCutOffs <- qnorm(1 - vProp)

  if (length(vCutOffs) == 1) {
    mMixRes[, binIDX] <- as.integer(mNormRes[, binIDX] >= vCutOffs)
  } else {
    mMixRes[, binIDX] <- t(apply(mNormRes[, vEPType == "Binary"], 1, function(x) {
      as.integer(x >= vCutOffs)
    }))
  }
  mMixRes
}
