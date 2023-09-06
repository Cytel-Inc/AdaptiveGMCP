
#Per look testing

PerLookMCPAnalysis<- function(mcpObj)
{
  P_Adj <- as.data.frame(apply(mcpObj$WH,1,Comb.Test,
                               p=mcpObj$p_raw ,cr=mcpObj$Correlation,upscale=T))
  names(P_Adj) <- paste('PAdj',mcpObj$CurrentLook, sep = '')

  PooledDF <- cbind(mcpObj$WH[, grep('H',names(mcpObj$WH))],P_Adj)

  if(mcpObj$CurrentLook == 1)
  {
    mcpObj$AdjPValues <- PooledDF
  }else
  {
    mcpObj$AdjPValues <- merge(mcpObj$AdjPValues, PooledDF, all = T)
  }

  if(mcpObj$CurrentLook > 1)
  {
    Comb_p <- c()
    for(i in 1:nrow(mcpObj$AdjPValues))
    {
      W_Norm <- mcpObj$W_Norm
      Comb_p[i] <- CombinedPvalue(CurrentLook = mcpObj$CurrentLook,
                                  adjPValue = mcpObj$AdjPValues[i,], W_Norm = W_Norm)
    }
    Comb_p <- as.data.frame(Comb_p)
    Comb_P_name <- paste('Comb_P',mcpObj$CurrentLook, sep = '')
    names(Comb_p) <- Comb_P_name
    mcpObj$AdjPValues <- cbind(mcpObj$AdjPValues, Comb_p)
  }


  for(idx in mcpObj$IndexSet) ## Closed test for each primary hypothesis
  {
    if(!mcpObj$rej_flag_Prev[idx]) #If not rejected in earlier look
    {

      if(mcpObj$CurrentLook == 1)
      {
        Intersect_IDX <- which(mcpObj$AdjPValues[idx]==1)
        mcpObj$rej_flag_Curr[idx] <- all(mcpObj$AdjPValues['PAdj1'][Intersect_IDX,] <= mcpObj$CutOff)

      }else #Combined P-Values
      {
        Intersect_IDX <- which(mcpObj$AdjPValues[idx]==1 &
                             !is.na(mcpObj$AdjPValues[paste('PAdj',mcpObj$CurrentLook,sep = '')]))

        mcpObj$rej_flag_Curr[idx] <- all( mcpObj$AdjPValues[Intersect_IDX,Comb_P_name] <= mcpObj$CutOff)

      }

      if(mcpObj$rej_flag_Curr[idx])
      {
        #If any hypothesis is rejected all the intersection hypothesis containing that can be removed
        mcpObj$WH <- mcpObj$WH[mcpObj$WH[idx] != 1,]
      }
    }
  }
  notRejected <- names(mcpObj$rej_flag_Curr[!mcpObj$rej_flag_Curr])
  Droped <- names(mcpObj$DropedFlag[mcpObj$DropedFlag])
  mcpObj$IndexSet <- setdiff(notRejected,Droped)
  mcpObj

}



#Combined P-value(Inverse Normal) assuming equal spacing

CombinedPvalue <- function(CurrentLook, adjPValue, W_Norm)
{
  weights <- rep(W_Norm[CurrentLook-1],CurrentLook)
  p_look <- as.numeric(adjPValue[,grep('PAdj',names(adjPValue))])
  if(any(is.na(p_look)))
  {
    return(NA)
  }else
  {
    1- pnorm(sum(weights*qnorm(1-p_look)))
  }

}


#Collect User Input Raw P-values
getRawPValues <- function(mcpObj)
{
  P_raw <- c()
  cat('User Input for the look : ',mcpObj$CurrentLook,'\n')
  for (i in mcpObj$IndexSet) {
    P_raw[i] <- as.numeric(readline(prompt = paste('Enter the raw P-Values for ', i, ' : ')))
  }
  P_raw
}

#User Input Selection
do_Selection <- function(mcpObj)
{
  SelectFlag <- readline(prompt = paste('Selection for the look ',(mcpObj$CurrentLook+1),' (y/n) :'))
  if(SelectFlag == 'y')
  {
    mcpObj$SelectionLook <- c(mcpObj$SelectionLook, (mcpObj$CurrentLook+1))

    cat('Out of ', paste(mcpObj$IndexSet, collapse = ', '),
        'Select the primary hypothesis to carry forward into the next look (e.g. H2,H3)','\n')

    select_index <- readline()

    mcpObj$SelectedIndex <- stringr::str_trim(
      unlist(strsplit(select_index,split = ',')),
      'both')

    for(i in 1:length(mcpObj$DropedFlag))
    {
      if(!mcpObj$DropedFlag[i] & !mcpObj$rej_flag_Prev[i])
        mcpObj$DropedFlag[i] = !any(mcpObj$SelectedIndex == names(mcpObj$DropedFlag[i]))
    }

    mcpObj$IndexSet <- intersect(mcpObj$IndexSet, mcpObj$SelectedIndex)
    mcpObj$WH <- mcpObj$WH[which(apply(mcpObj$WH[mcpObj$IndexSet], 1, sum, na.rm=T) != 0),]
  }
  mcpObj
}


get_numeric_part <- function(vec) {
  numeric_part <- as.numeric(gsub("[^0-9]", "", vec))
  return(numeric_part)
}

#Modification of weights and graph for the continuing hypothesis
do_modifyStrategy <- function(mcpObj, showExistingStrategy = T)
{
  ModificationFlag <- readline(prompt = paste('Change the testing Strategy from the look :',(mcpObj$CurrentLook+1),' (y/n) : \n'))

  if(ModificationFlag == 'y')
  {
    if(showExistingStrategy)
    {
      cat('Existing Strategy for reference \n')
      print(mcpObj$WH)
    }

    mcpObj$ModificationLook <- c(mcpObj$ModificationLook,mcpObj$CurrentLook)

    #User input weights
    cat("Enter the new weights for (",paste(mcpObj$IndexSet, collapse = ', '),') as comma seperated values (e.g.- 0.5,0.5) :\n')
    new_weights <- readline()

    mcpObj$newWeights <- as.numeric(stringr::str_trim(
      unlist(strsplit(new_weights,split = ',')),
      'both'))
    names(mcpObj$newWeights) = paste('Weight',get_numeric_part(mcpObj$IndexSet),sep='')

    #User input transition matrix
    m = length(mcpObj$IndexSet)
    new_G = matrix(0, nrow = m, ncol = m)
    cat("\n Enter the elements of the transition matrix G=(gij) \n")
    for(i in 1:m)
    {
      for (j in 1:m) {
        if(i != j)
        {
          cat('Enter g(',paste(mcpObj$IndexSet[i],mcpObj$IndexSet[j],sep = ','),') :\n')
          new_G[i,j] = as.numeric(readline())
        }
      }
    }
    colnames(new_G) <- mcpObj$IndexSet
    mcpObj$newG <- new_G

    modifiedWeights <- modifyIntersectWeights(mcpObj)
    return(modifiedWeights$mcpObj)
  }else
  {
    return(mcpObj)
  }
}

#Merge Existing and Modified weights
modifyIntersectWeights <- function(mcpObj)
{
  A <- mcpObj$WH
  B <- as.data.frame(generateWeights(g = mcpObj$newG,
                                     w = mcpObj$newWeights))
  names(B) <- c(mcpObj$IndexSet,
                          paste('Weight',get_numeric_part(mcpObj$IndexSet),sep=''))

  ## If the weights for the hypothesis belongs to IA set(Section 8.2) needs to be
  ## considered from the before modification strategy

  isIAsameIC = F # Ajoy.M: it should consider the earlier strategy(4thSept,23)

  if(isIAsameIC == F)
  {
    HypNotAvail <- setdiff(names(A)[ grep('H',names(A))], mcpObj$IndexSet)
    WNotAvail <- paste('Weight',get_numeric_part(HypNotAvail),sep='')
    if(length(HypNotAvail) < 2 )
    {
      B[HypNotAvail] <- rep(0, nrow(B))
      B[WNotAvail] <- rep(0, nrow(B))
    }else
    {
      B[HypNotAvail] <- matrix(0,nrow = nrow(B), ncol = length(HypNotAvail))
      B[WNotAvail] <- matrix(0,nrow = nrow(B), ncol = length(WNotAvail))
    }
    UpdatedWeights <- A
    B <- B[,colnames(A)] #Reordering columns of B
    HA <- apply(A[,grep('H',names(A))], 1, paste,collapse = '')
    HB <- apply(B[,grep('H',names(B))], 1, paste,collapse = '')

    if(!all(HB %in% HA)) return("Error in MergeWeights")
    modrows <- c()

    for (i in 1:nrow(B)) {
      r <- which(HA==HB[i])
      modrows <- c(modrows, r)

      l1 <- A[r, grep('Weight',names(A))]
      l2 <- B[i, grep('Weight',names(B))]
      for(l in names(l1))
      {
        if(l %in% names(l2))
        {
          UpdatedWeights[r,l] = l2[l]
        }else
        {
          UpdatedWeights[r,l] = NA
        }
      }
    }

    mcpObj$WH <- UpdatedWeights
    return(list('mcpObj'=mcpObj, 'modifiedColumns'=modrows))
  }else
  {
    mcpObj$WH <- B
    return(list('mcpObj'=mcpObj))
  }

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
  cat('Analysis results for Look : ',mcpObj$CurrentLook,'\n')
  cat('\n')
  cat('\n')

  cat('Weights for the intersection hypothesis at Look : ',mcpObj$CurrentLook,'\n')
  print(mcpObj$WH_Prev)

  cat('\n')
  cat('\n')

  cat('Adj P-values for the intersection hypothesis at Look : ',mcpObj$CurrentLook,'\n')
  print(mcpObj$AdjPValues)

  cat('\n')
  cat('\n')

  cat('Rejection Status of primary Hypothesis at Look : ',mcpObj$CurrentLook,'\n')
  print(mcpObj$rej_flag_Curr)

  cat('\n')
  cat('\n')
}

#Stop Trial
StopTrial <- function(mcpObj)
{
  if(all(mcpObj$rej_flag_Curr[!mcpObj$DropedFlag]==T))#Stop if all the primary hypothesis are rejected
  {
    return(T)

  }else if(mcpObj$CurrentLook == mcpObj$LastLook)
  {
    return(T)
  }else
  {
    F
  }
}
