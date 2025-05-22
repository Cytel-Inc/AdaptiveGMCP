# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

#Combination P-value Method
SingleSimCombPValue2 <- function(simID, gmcpSimObj, preSimObjs){
  tryCatch(
    SingleSimCombPValue(simID, gmcpSimObj, preSimObjs),
    error = function(e){
      sprintf("Error Simulation %d ", simID)
    }
  )
}

#CER Method
SingleSimCER2 <- function(simID, gmcpSimObj, preSimObjs) {
  tryCatch(
    SingleSimCER(simID, gmcpSimObj, preSimObjs),
    error = function(e) {
      sprintf("Error in Simulation %d: %s", simID, e$message)
    }
  )
}

