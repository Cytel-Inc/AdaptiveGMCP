# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

########### Identification of disjoint Sets from the intersection Hypothesis ###########
# J : binary index for the intersection hypothesis
# test.type : planned test type
# HypoMap : dataframe to identify the endpoint specific hypothesis
#------------------------------------------------------- -
# Connected Subsets
connSets <- function(J, w, test.type, HypoMap) {
  if (test.type == "Partly-Parametric" || test.type == "Parametric") {
    df1 <- HypoMap[which(J == 1), ]
    unique_elements <- unique(df1$Groups)
    disSets <- lapply(
      unique_elements,
      function(x) which(HypoMap$Groups == x & J != 0 & w != 0)
    )
    disSets <- disSets[lapply(disSets, length) > 0]

    lenSets <- unlist(lapply(disSets, length))
    if (length(lenSets) == 1 & all(lenSets > 1)) {
      Method <- "Parametric"
    } else if (length(lenSets) > 1 & any(lenSets > 1)) {
      Method <- "Mixed"
    } else {
      Method <- "Non_Parametric"
    }
    return(list("connSets" = disSets, "Method" = Method))
  } else {
    hyp <- 1:length(J)
    disSets <- lapply(hyp[J == 1], function(x) x)
    return(list("connSets" = disSets, "Method" = "Non_Parametric"))
  }
}
#------------------------------------------------------------- -

# Connected parametric subsets
connParamNParmSets <- function(conn_Sets) {
  lSets <- lapply(conn_Sets, length)
  Param_grp <- conn_Sets[which(lSets > 1)]
  nParam_grp <- unlist(conn_Sets[which(lSets == 1)])
  list("ParamGrps" = Param_grp, "NParamGrp" = nParam_grp)
}


#------------------------------------------------------------- -
# Create mapping between hypothesis arms and endpoints
getHypoMap <- function(des.type, nHypothesis, nEps, nArms) {
  if (des.type == "MAMSMEP") {
    HypoMap <- data.frame(
      "Hypothesis" = paste("H", 1:nHypothesis, sep = ""),
      "Groups" = rep(1:nEps, each = (nArms - 1)),
      "Control" = rep(1, nHypothesis),
      "Treatment" = rep(2:nArms, nEps),
      row.names = NULL
    )
  }
  HypoMap
}

#------------------------------------------------------------- -
# Create mapping between hypothesis arms and endpoints
getHypoMap2 <- function(des.type, nHypothesis, nEps, nArms, lEpType) {
  if (des.type == "MAMSMEP") {
    HypoMap <- data.frame(
      "Hypothesis" = paste("H", 1:nHypothesis, sep = ""),
      "Groups" = rep(1:nEps, each = (nArms - 1)),
      "EpType" = rep(unlist(lEpType), each = (nArms - 1)),
      "Control" = rep(1, nHypothesis),
      "Treatment" = rep(2:nArms, nEps),
      row.names = NULL
    )
  }
  HypoMap
}



#-------------- -
# Get Arms from Hypothesis index
#-------------- -
getArmsFromHypo <- function(SetH, Hypo_map) {
  if (length(SetH) == 0) {
    return(NULL)
  } else {
    df1 <- Hypo_map[Hypo_map$Hypo %in% SetH, ]
    return(unique(c(df1$Control, df1$Treatment)))
  }
}
