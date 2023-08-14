

#' Function to perform Adaptive GMCP by using Combined P-Value(Inverse Normal) by user input P-values
#' @param K  Number of looks(default = 3)
#' @param D  Number of Hypothesis(default = 4)
#' @param WI Vector of Initial Weights for Global Null(default = equal weights)
#' @param G  Transition Matrix (default = equal distribution)
#' @param Threshold Vector of Critical points to reject H0i in terms of p-values in the interim looks
#' @param Correlation Matrix of correlation between test statistics, NA if the correlation is unknown
#' @param ShowGraph Logical, If set as True graphs will be printed


adaptGMCP_CP <- function(
    K = 3, D = 4,
    WI = rep(1/4,4),
    G = matrix(c(0,1/3,1/3,1/3,  1/3,0,1/3,1/3, 1/3,1/3,0,1/3, 1/3,1/3,1/3,0),nrow = 4),
    Threshold = c(0.001,0.0015,0.025),
    Correlation = matrix(c(1,0.5,NA,NA,  0.5,1,NA,NA,  NA,NA,1,NA, NA,NA,NA,1),nrow = 4),
    ShowGraph = F
)
{
  GlobalIndexSet <- paste("H",1:D, sep ='')

  ## Weights for all intersection hypothesis
  WH <- as.data.frame(generateWeights(G, WI))
  colnames(WH) <- c(GlobalIndexSet, paste("Weight",1:D, sep = ''))

  rej_flag_Prev <- rej_flag_Curr <-  rep(FALSE, D)
  names(rej_flag_Prev) <- names(rej_flag_Curr) <- paste("H",1:D, sep = '')


  mcpObj <- list('CurrentLook'= NA, 'IndexSet'= GlobalIndexSet,
                 'p_raw' = NA, 'WH'= WH,
                 'Correlation' = Correlation,'AdjPValues' = NA,
                 'CutOff' = NA, 'rej_flag_Prev'= rej_flag_Prev , 'rej_flag_Curr'= rej_flag_Curr)

  SummaryStatLs <- list('LookID'=NA, 'Rej_Flag'=NA)
  for(j in 1:K) ## Loop for each Interim Analysis
  {
    mcpObj$CurrentLook <- j
    p_raw <- getRawPValues(mcpObj) ## User Input raw p-values

    mcpObj$p_raw <- addNAPvalue(p_raw, GlobalIndexSet)
    mcpObj$CutOff <- Threshold[j]

    mcpObj <- PerLookMCPAnalysis(mcpObj)
    mcpObj$rej_flag_Prev <- mcpObj$rej_flag_Curr

    ShowResults(mcpObj)
    if(all(mcpObj$rej_flag_Curr==T)) break
  }
}




PerLookMCPAnalysis<- function(mcpObj)
{
  p_intersect <- t(apply(mcpObj$WH,1,pvals.dunnett,
                         p=mcpObj$p_raw ,cr=mcpObj$Correlation,upscale=T))

  P_Adj <- as.data.frame(apply(p_intersect, 1, min, na.rm=T))
  names(P_Adj) <- paste('PAdj',mcpObj$CurrentLook, sep = '')

  PooledDF <- cbind(mcpObj$WH[, grep('H',names(mcpObj$WH))],P_Adj)

  if(mcpObj$CurrentLook == 1)
  {
    mcpObj$AdjPValues <- PooledDF
  }else
  {
    mcpObj$AdjPValues <- merge(mcpObj$AdjPValues, PooledDF, all = T)
  }

  for(idx in mcpObj$IndexSet) ## Closed test for each primary hypothesis
  {
    if(!mcpObj$rej_flag_Prev[idx]) #If not rejected in earlier look
    {

      if(mcpObj$CurrentLook == 1)
      {
        Intersect_IDX <- which(mcpObj$AdjPValues[idx]==1)
        mcpObj$rej_flag_Curr[idx] <- all(mcpObj$AdjPValues['PAdj1'][Intersect_IDX,] < mcpObj$CutOff)

      }else #Combined P-Values
      {
        Intersect_IDX <- which(mcpObj$AdjPValues[idx]==1 &
                                 !is.na(mcpObj$AdjPValues[paste('PAdj',mcpObj$CurrentLook,sep = '')]))
        ComPVal <- c()
        for(h in Intersect_IDX)
        {
          ComPVal <- c(ComPVal,
                       CombinedPvalue(CurrentLook = mcpObj$CurrentLook,
                                      adjPValue = mcpObj$AdjPValues[h,]))
        }

        mcpObj$rej_flag_Curr[idx] <- all(ComPVal < mcpObj$CutOff)

      }

      if(mcpObj$rej_flag_Curr[idx])
      {
        #If any hypothesis is rejected all the intersection hypothesis containing that can be removed
        mcpObj$WH <- mcpObj$WH[mcpObj$WH[idx] != 1,]
      }
    }
  }
  notRejected <- names(mcpObj$rej_flag_Curr[!mcpObj$rej_flag_Curr])
  mcpObj$IndexSet <- notRejected

  mcpObj

}

#Combined P-value(Inverse Normal) assuming equal spacing

CombinedPvalue <- function(CurrentLook, adjPValue)
{
  weights <- rep(sqrt(1/CurrentLook),CurrentLook)
  p_look <- as.numeric(adjPValue[,grep('PAdj',names(adjPValue))])
  1- pnorm(sum(weights*qnorm(1-p_look)))
}


#Collect User Input Raw P-values
getRawPValues <- function(mcpObj)
{
  P_raw <- c()
  cat('User Input for Look : ',mcpObj$CurrentLook,'\n')
  for (i in mcpObj$IndexSet) {
    P_raw[i] <- as.numeric(readline(prompt = paste('Enter Raw P-Values for ', i, ' : ')))
  }
  P_raw
}

#Add NA for already rejected hypothesis
addNAPvalue <- function(p_raw, GlobalIndexSet)
{
  p_raw_na <- rep(NA,length(GlobalIndexSet))
  names(p_raw_na) <- GlobalIndexSet
  p_raw_na[names(p_raw)] <- p_raw
  p_raw_na
}


#Console Output of Look Wise results
ShowResults <- function(mcpObj)
{
  cat('\n')
  cat('\n')

  cat('Intersection hypothesis at Look : ',mcpObj$CurrentLook,'\n')
  print(mcpObj$AdjPValues)

  cat('\n')
  cat('\n')

  cat('Rejection Status of primary Hypothesis at Look : ',mcpObj$CurrentLook,'\n')
  print(mcpObj$rej_flag_Curr)

  cat('\n')
  cat('\n')
}










