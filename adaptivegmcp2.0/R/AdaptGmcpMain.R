

#' Function to perform Adaptive GMCP by using Combined P-Value(Inverse Normal) by user input P-values
#' @param K  Number of looks(default = 3)
#' @param D  Number of Hypothesis(default = 4)
#' @param WI Vector of Initial Weights for Global Null(default = rep(1/4,4))
#' @param G  Transition Matrix (default = matrix(c(0,1/3,1/3,1/3,  1/3,0,1/3,1/3, 1/3,1/3,0,1/3, 1/3,1/3,1/3,0),nrow = 4))
#' @param Threshold Vector of Critical points to reject H0i in terms of p-values in the interim looks (default = c(0.001,0.0015,0.025))
#' @param Correlation Matrix of correlation between test statistics, NA if the correlation is unknown
#' @param ShowGraph Logical, If set as True graphs will be printed


adaptGMCP_PC <- function(
    K = 3, D = 4,
    WI = rep(1/4,4),
    G = matrix(c(0,1/3,1/3,1/3,  1/3,0,1/3,1/3, 1/3,1/3,0,1/3, 1/3,1/3,1/3,0),nrow = 4),
    Threshold = c(0.001,0.0015,0.025),
    Correlation = matrix(c(1,0.5,NA,NA,  0.5,1,NA,NA,  NA,NA,1,NA, NA,NA,NA,1),nrow = 4),
    Selection=T,
    UpdateStrategy=T,
    ShowGraph = F
)
{
  GlobalIndexSet <- paste("H",1:D, sep ='')

  ## Weights for all intersection hypothesis
  WH <- as.data.frame(generateWeights(G, WI))
  colnames(WH) <- c(GlobalIndexSet, paste("Weight",1:D, sep = ''))

  rej_flag_Prev <- rej_flag_Curr <- DropedFlag <- rep(FALSE, D)
  names(rej_flag_Prev) <- names(rej_flag_Curr) <- names(DropedFlag) <-paste("H",1:D, sep = '')

  #Modified Weights and Graph based on interim data
  modifiedStrategy <- list('Modify'=F,
                           'ModificationLook'=c(),
                           'newWeights' =NA,
                           'newG' = NA)

  mcpObj <- list('CurrentLook'= NA,
                 'IntialWeights'=WI,
                 'IntialHypothesis'=GlobalIndexSet,
                 'IndexSet'= GlobalIndexSet,
                 'p_raw' = NA,
                 'WH'= WH,
                 'Correlation' = Correlation,
                 'AdjPValues' = NA,
                 'CutOff' = NA,
                 'rej_flag_Prev'= rej_flag_Prev ,
                 'rej_flag_Curr'= rej_flag_Curr,
                 'SelectionLook'= c(),
                 'SelectedIndex'=NA,
                 'DropedFlag'=DropedFlag,
                 'LastLook'=K,
                 'ModifiedStrategy'= modifiedStrategy)


  look = 1; ContTrial = T

  while(ContTrial) ## Loop for each Interim Analysis
  {
    mcpObj$CurrentLook <- look
    p_raw <- getRawPValues(mcpObj) ## User Input raw p-values

    mcpObj$p_raw <- addNAPvalue(p_raw, GlobalIndexSet)
    mcpObj$CutOff <- Threshold[look]

    mcpObj <- PerLookMCPAnalysis(mcpObj)
    mcpObj$rej_flag_Prev <- mcpObj$rej_flag_Curr

    ShowResults(mcpObj)

    if(StopTrial(mcpObj)) break  #Early Stoping

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

    look <- look + 1
  }

  #mcpObj
}








