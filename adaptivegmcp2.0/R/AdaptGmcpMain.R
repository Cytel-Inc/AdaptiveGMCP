

#' Function to perform Adaptive GMCP by using Combined P-Value(Inverse Normal) by user input P-values
#' @param WI Vector of Initial Weights for Global Null(default = \code{rep(1/4,4)})
#' @param G  Transition Matrix (default = \code{matrix(c(0,1/3,1/3,1/3,  1/3,0,1/3,1/3, 1/3,1/3,0,1/3, 1/3,1/3,1/3,0),nrow = 4)})
#' @param Threshold Vector of Critical points to reject H0i in terms of p-values in the interim looks (default = \code{c(0.001,0.0015,0.025)})
#' @param Correlation Matrix of correlation between test statistics, NA if the correlation is unknown
#' @param W_Norm Weights for inverse normal combining p-value function(default = \code{sqrt(c(1/K, 2/K,..,(K-1)/K)})
#' @Selection Logical: True if selection required at interim(default = False)
#' @UpdateStrategy Logical: True if modification of weights and testing strategy is required at interim(default = False)
#' @export
adaptGMCP_PC <- function(
    WI = rep(1/4,4),
    G = matrix(c(0,1/3,1/3,1/3,  1/3,0,1/3,1/3, 1/3,1/3,0,1/3, 1/3,1/3,1/3,0),nrow = 4),
    Threshold = c(0.001,0.0015,0.025),
    Correlation = matrix(c(1,0.5,NA,NA,  0.5,1,NA,NA,  NA,NA,1,NA, NA,NA,NA,1),nrow = 4),
    W_Norm = sqrt(1:(length(Threshold)-1)/length(Threshold)),
    Selection=F,
    UpdateStrategy=F
)
{
  D <- length(WI)
  K <- length(Threshold)
  GlobalIndexSet <- paste("H",1:D, sep ='')

  ## Weights for all intersection hypothesis
  WH <- as.data.frame(generateWeights(G, WI))
  colnames(WH) <- c(GlobalIndexSet, paste("Weight",1:D, sep = ''))

  rej_flag_Prev <- rej_flag_Curr <- DropedFlag <- rep(FALSE, D)
  names(rej_flag_Prev) <- names(rej_flag_Curr) <- names(DropedFlag) <-paste("H",1:D, sep = '')


  #info to run per look analysis
  mcpObj <- list('CurrentLook'= NA,
                 'IntialWeights'=WI,
                 'IntialHypothesis'=GlobalIndexSet,
                 'IndexSet'= GlobalIndexSet,
                 'p_raw' = NA,
                 'WH_Prev' = WH,
                 'WH'= WH,
                 'Correlation' = Correlation,
                 'AdjPValues' = NA,
                 'W_Norm' = W_Norm,
                 'CutOff' = NA,
                 'rej_flag_Prev'= rej_flag_Prev ,
                 'rej_flag_Curr'= rej_flag_Curr,
                 'SelectionLook'= c(),
                 'SelectedIndex'=NA,
                 'DropedFlag'=DropedFlag,
                 'LastLook'=K,
                 'Modify'=F,
                 'ModificationLook'=c(),
                 'newWeights' =NA,
                 'newG' = NA
                 )

  #to Store all look info
  allInfo <- list()

  look = 1; ContTrial = T

  while(ContTrial) ## Loop for each Interim Analysis
  {
    mcpObj$CurrentLook <- look
    p_raw <- getRawPValues(mcpObj) ## User Input raw p-values

    mcpObj$p_raw <- addNAPvalue(p_raw, GlobalIndexSet)
    mcpObj$CutOff <- Threshold[look]

    mcpObj <- PerLookMCPAnalysis(mcpObj)
    ShowResults(mcpObj) #Current look information

    lookID <- paste('Look',look,sep = '')
    allInfo[[lookID]] <- mcpObj

    mcpObj$rej_flag_Prev <- mcpObj$rej_flag_Curr
    mcpObj$WH_Prev <- mcpObj$WH

    if(StopTrial(mcpObj))
    {
      break  #Early Stoping

    }else  #Next look
    {
      # Selection for next look
      if(Selection & (look< K))
      {
        mcpObj <- do_Selection(mcpObj)
      }

      #Modify the weights and testing strategy
      if(UpdateStrategy & (look< K) & (length(mcpObj$IndexSet)>1))
      {
        mcpObj <- do_modifyStrategy(mcpObj)
      }
    }
    look <- look + 1
  }
  #allInfo
}








