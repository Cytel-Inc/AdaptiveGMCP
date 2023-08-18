
#Combined P-value(Inverse Normal) assuming equal spacing

CombinedPvalue <- function(CurrentLook, adjPValue)
{
  weights <- rep(sqrt(1/CurrentLook),CurrentLook)
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
  cat('User Input for Look : ',mcpObj$CurrentLook,'\n')
  for (i in mcpObj$IndexSet) {
    P_raw[i] <- as.numeric(readline(prompt = paste('Enter Raw P-Values for ', i, ' : ')))
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
