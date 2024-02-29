

#Perform Single Simulation for combining p-values method
#' @param simID Simulation ID
#' @param gmcpSimObj obj with simulation inputs
#' @param preSimObjs obj with intermediate inputs
SingleSimCombPValue <- function(simID, gmcpSimObj,preSimObjs)
{
  #Initialize Intermediate Inputs
  mcpObj <- initialize_mcpObj(gmcpSimObj = gmcpSimObj, preSimObjs = preSimObjs)
  SummStatDF <- mcpObj$SummStatBlank
  ArmWiseDF <- mcpObj$ArmWiseDataBlank

  #look-wise test
  while(mcpObj$ContTrial)
  {
    #Get the per arm incremental sample size
    Arms.SS.Incr <- getInterimSSIncr(lookID=mcpObj$CurrentLook,
                                     PlanSSIncr=mcpObj$planSS$IncrementalSamples,
                                     ArmsPresent=mcpObj$ArmsPresent,
                                     ArmsRetained=mcpObj$ArmsRetained,
                                     Arms.alloc.ratio=mcpObj$Arms.alloc.ratio,
                                     ImplicitSSR=mcpObj$ImplicitSSR)

    #Get incremental summary
    currLookDataIncr <- genIncrLookSummaryDOM(SimSeed=mcpObj$SimSeed,
                                          simID=simID,
                                          lookID=mcpObj$CurrentLook,
                                          Arms.Mean=mcpObj$Arms.Mean,
                                          Arms.std.dev=mcpObj$Arms.std.dev,
                                          Arms.alloc.ratio=mcpObj$Arms.alloc.ratio,
                                          Arms.SS=Arms.SS.Incr,
                                          EPCorr=mcpObj$EP.Corr,
                                          ArmsPresent=mcpObj$ArmsPresent,
                                          HypoPresent=mcpObj$HypoPresent,
                                          HypoMap=mcpObj$HypoMap)
    ArmData <- currLookDataIncr$ArmData

    #Compute Test Stat(Incr.) and p-values(Incr.)
    SummStat <- getPerLookTestStatDOM(simID = simID,
                                      lookID=mcpObj$CurrentLook,
                                      TestStat = mcpObj$TestStat,
                                      Arms.std.dev=mcpObj$Arms.std.dev,
                                      IncrLookSummaryDOM=currLookDataIncr,
                                      HypoMap=mcpObj$HypoMap,
                                      Cumulative=FALSE)
    #Perform per look Test
    mcpObj <- perLookTest(Arms.SS.Incr=Arms.SS.Incr,SummStat = SummStat, mcpObj=mcpObj)

    #Available Arms & hypothesis after rejection
    mcpObj$ArmsPresent <- getArmsPresent(ArmsPresent = mcpObj$ArmsPresent, rejflags = mcpObj$rej_flag_Curr, HypoMap = mcpObj$HypoMap)
    mcpObj$HypoPresent <- !mcpObj$rej_flag_Curr
    mcpObj$IndexSet <- mcpObj$HypoMap$Hypothesis[mcpObj$HypoPresent]

    rejStatus <- data.frame(matrix(mcpObj$rej_flag_Curr, nrow = 1))
    names(rejStatus) <- paste('RejStatus',get_numeric_part(mcpObj$HypoMap$Hypothesis),sep='')

    SummStat <- plyr::rbind.fill(mcpObj$SummStatBlank,SummStat)
    SummStat <- fillNa(1,SummStat,rejStatus)
    SummStatDF <- plyr::rbind.fill(SummStatDF,SummStat)
    ArmWiseDF <- plyr::rbind.fill(ArmWiseDF, ArmData)
    mcpObj$SummStatDF <- SummStatDF
    mcpObj$ArmDataDF <- ArmWiseDF

    mcpObj$rej_flag_Prev <- mcpObj$rej_flag_Curr
    mcpObj$WH_Prev <- mcpObj$WH

    #Check for Early Stopping Conditions
    if(StopTrial(mcpObj))
    {
      mcpObj$ContTrial = F
      break  #Early Stopping

    }else  #Next look
    {
      # Selection for next look
      if(mcpObj$Selection & (mcpObj$CurrentLook < mcpObj$LastLook))
      {
        mcpObj <- do_SelectionSim2(mcpObj)
      }
    }
    mcpObj$CurrentLook  <- mcpObj$CurrentLook + 1

    # End of Look-wise While Loop #
  }

  #Power Table
  powerCountDF <- CountPower(simID = simID, SummaryStatFile = SummStatDF, TrueNull = mcpObj$TrueNull)
  EffCountDF <- CountEfficacy(simID = simID, SummaryStatFile = SummStatDF)

  SelectionDF <- data.frame()
  #Selection Summary Table
  if(mcpObj$Selection & length(mcpObj$SelectedIndex)>0)
  {
    SelectionDF <-data.frame('simID'=simID,
                             'SelectedHypothesis'=mcpObj$SelectedIndex)
  }

  list('SummStatDF'=SummStatDF,
       'ArmWiseDF'=ArmWiseDF,
       'powerCountDF'=powerCountDF,
       'EfficacyTable'=EffCountDF,
       'SelectionDF'=SelectionDF)
}


#Perform Single Simulation for CER method
#' @param simID Simulation ID
#' @param gmcpSimObj obj with simulation inputs
#' @param preSimObjs obj with intermediate inputs
SingleSimCER <- function(simID, gmcpSimObj,preSimObjs)
{
  #Initialize Intermediate Inputs
  mcpObj <- initialize_mcpObj(gmcpSimObj = gmcpSimObj, preSimObjs = preSimObjs)
  SummStatDF <- mcpObj$SummStatBlank
  ArmWiseDF <- mcpObj$ArmWiseDataBlank

  #look-wise test
  while(mcpObj$ContTrial)
  {
    #Get the per arm incremental sample size
    Arms.SS.Incr <- getInterimSSIncr(lookID=mcpObj$CurrentLook,
                                     PlanSSIncr=mcpObj$planSS$IncrementalSamples,
                                     ArmsPresent=mcpObj$ArmsPresent,
                                     ArmsRetained=mcpObj$ArmsRetained,
                                     Arms.alloc.ratio=mcpObj$Arms.alloc.ratio,
                                     ImplicitSSR=mcpObj$ImplicitSSR)

    #If current Stage = 2, check for adaptations
    if(mcpObj$CurrentLook==2)
    {
      PlanSSLk <- mcpObj$planSS$IncrementalSamples[mcpObj$CurrentLook,]
      CurrSSLk <- Arms.SS.Incr
      mcpObj$AdaptStage2 <- F

      #If Current look Sample Size is different than the plan Sample size then Adapt
      if(!all(PlanSSLk == CurrSSLk,na.rm = T))
      {
        mcpObj$AdaptStage2 <- T
        Stage2AllocSampleSize <- unlist(mcpObj$planSS$IncrementalSamples[1,])+
          unlist(CurrSSLk)
        mcpObj$Stage2AllocSampleSize <- rbind(mcpObj$planSS$IncrementalSamples[1,],
                                              Stage2AllocSampleSize)
        mcpObj$Stage2allocRatio <- unlist(CurrSSLk)/unlist(CurrSSLk[1])

        #Compute the Stage2 adaptive boundaries
        AdaptResults <- adaptBdryCER(mcpObj)
        mcpObj$AdaptObj <- AdaptResults
      }
    }


    #Generate look-wise incremental data
    currLookDataIncr <- genIncrLookSummaryDOM(SimSeed=mcpObj$SimSeed,
                                              simID=simID,
                                              lookID=mcpObj$CurrentLook,
                                              Arms.Mean=mcpObj$Arms.Mean,
                                              Arms.std.dev=mcpObj$Arms.std.dev,
                                              Arms.alloc.ratio=mcpObj$Arms.alloc.ratio,
                                              Arms.SS=Arms.SS.Incr,
                                              EPCorr=mcpObj$EP.Corr,
                                              ArmsPresent=mcpObj$ArmsPresent,
                                              HypoPresent=mcpObj$HypoPresent,
                                              HypoMap=mcpObj$HypoMap)
    ArmData <- currLookDataIncr$ArmData


    if(mcpObj$CurrentLook==1)
    {
      # Summary Statistics based on first look data
      SummStat <- getPerLookTestStatDOM(simID = simID,
                                        lookID = mcpObj$CurrentLook,
                                        TestStat = mcpObj$TestStat,
                                        Arms.std.dev=mcpObj$Arms.std.dev,
                                        IncrLookSummaryDOM = currLookDataIncr,
                                        HypoMap=mcpObj$HypoMap,
                                        Cumulative=FALSE)

    }else
    {
      if(mcpObj$FWERControl == 'CombinationTest')
      {
        # Summary Statistics based on 2nd look incremental data
        SummStat <- getPerLookTestStatDOM(simID = simID,
                                          lookID=mcpObj$CurrentLook,
                                          TestStat = mcpObj$TestStat,
                                          Arms.std.dev=mcpObj$Arms.std.dev,
                                          IncrLookSummaryDOM=currLookDataIncr,
                                          HypoMap=mcpObj$HypoMap,
                                          Cumulative=FALSE)

        # Stage-1 raw p-values(Incr.)
        pValIncrPrev <- mcpObj$SummStatDF[mcpObj$SummStatDF$LookID == (mcpObj$CurrentLook-1),
                                          grep('RawPvalues', names(mcpObj$SummStatDF))]
        # Stage-2 raw p-values(Incr.)
        pValIncrCurr <- SummStat[,grep('RawPvalues', names(SummStat))]

        pValIncr <- plyr::rbind.fill(pValIncrPrev,pValIncrCurr)

        # Converted to Z Statistics(Incr.)
        zIncr <- apply(pValIncr, 2, function(x){qnorm(1-x)})

        W_Norm <- mcpObj$InvNormWeights$W_Norm
        if(is.vector(W_Norm) & mcpObj$CurrentLook == 2)
        {
          #Inverse Normal Weights for two looks based on pre-planed look positions
          W_Inv <- W_Norm
        }else if(is.matrix(W_Norm)) #W_Norm is a matrix for more than 2 looks
        {
          #Inverse Normal Weights for more than two looks  based on pre-planed look positions
          W_Inv <- W_Norm[((mcpObj$CurrentLook)-1),1:(mcpObj$CurrentLook)]
        }else
        {
          #Error in CombinedPvalue function
          return('Error in CombinedPvalue function')
        }

        #Test The computed Inverse Normal Weights
        if(abs(sum(W_Inv^2)-1) > 1E-6) stop('Error: abs(sum(W_Inv^2)-1) < 1E-6 not true| function: CombinedPvalue')

        #Inverse Normal Combination
        combStage2TestStat <- unlist(lapply(1:ncol(zIncr), function(i){sum(W_Inv*zIncr[,i])}))

        combStage2pVal <- unlist(lapply(1:ncol(zIncr), function(i){
          1-pnorm(sum(W_Inv*zIncr[,i]))
        }))
        SummStat[,grep('TestStat', names(SummStat))] <- combStage2TestStat
        SummStat[,grep('RawPvalues', names(SummStat))] <- combStage2pVal

      }else
      {
        # Summary Statistics based on 2nd look cumulative data
        SummStat <- getPerLookTestStatDOM(simID = simID,
                                          lookID=mcpObj$CurrentLook,
                                          TestStat = mcpObj$TestStat,
                                          Arms.std.dev=mcpObj$Arms.std.dev,
                                          IncrLookSummaryDOM=currLookDataIncr,
                                          IncrLookSummaryDOMPrev = IncrLookSummaryDOMPrev,
                                          HypoMap=mcpObj$HypoMap,
                                          Cumulative=TRUE)
      }

    }
    #Storing the current look incremental data to compute next look cumulative data
    IncrLookSummaryDOMPrev <- currLookDataIncr

    #Perform per look Test
    mcpObj <- perLookTest(Arms.SS.Incr=Arms.SS.Incr,SummStat = SummStat, mcpObj=mcpObj)

    #Available Arms & hypothesis after rejection
    mcpObj$ArmsPresent <- getArmsPresent(ArmsPresent = mcpObj$ArmsPresent, rejflags = mcpObj$rej_flag_Curr, HypoMap = mcpObj$HypoMap)
    mcpObj$HypoPresent <- !mcpObj$rej_flag_Curr
    mcpObj$IndexSet <- mcpObj$HypoMap$Hypothesis[mcpObj$HypoPresent]

    rejStatus <- data.frame(matrix(mcpObj$rej_flag_Curr, nrow = 1))
    names(rejStatus) <- paste('RejStatus',get_numeric_part(mcpObj$HypoMap$Hypothesis),sep='')

    SummStat <- plyr::rbind.fill(mcpObj$SummStatBlank,SummStat)
    SummStat <- fillNa(1,SummStat,rejStatus)
    SummStatDF <- plyr::rbind.fill(SummStatDF,SummStat)
    ArmWiseDF <- plyr::rbind.fill(ArmWiseDF, ArmData)
    mcpObj$SummStatDF <- SummStatDF
    mcpObj$ArmDataDF <- ArmWiseDF

    mcpObj$rej_flag_Prev <- mcpObj$rej_flag_Curr


    #Check for Early Stopping Conditions
    if(StopTrial(mcpObj))
    {
      mcpObj$ContTrial = F
      break  #Stop the Trial

    }else  #Next look
    {
      # Selection for next look
      if(mcpObj$Selection & (mcpObj$CurrentLook < mcpObj$LastLook))
      {
        mcpObj <- do_SelectionSim2(mcpObj)
      }

    }
    # move to next look
    mcpObj$CurrentLook  <- mcpObj$CurrentLook + 1

    # End of Look-wise Loop #
  }
  #Power Table
  powerCountDF <- CountPower(simID = simID, SummaryStatFile = SummStatDF, TrueNull = mcpObj$TrueNull)
  EffCountDF <- CountEfficacy(simID = simID, SummaryStatFile = SummStatDF)

  SelectionDF <- data.frame()
  #Selection Summary Table
  if(mcpObj$Selection & length(mcpObj$SelectedIndex)>0)
  {
    SelectionDF <-data.frame('simID'=simID,
                             'SelectedHypothesis'=mcpObj$SelectedIndex)
  }

  list('SummStatDF'=SummStatDF,
       'ArmWiseDF'=ArmWiseDF,
       'powerCountDF'=powerCountDF,
       'EfficacyTable'=EffCountDF,
       'SelectionDF'=SelectionDF)
}






#Initialization of mcpObj to run look-wise analysis
initialize_mcpObj <- function(gmcpSimObj, preSimObjs)
{
  SummStatNames <- c('SimID','LookID',
                     paste('Delta',1:gmcpSimObj$nHypothesis,sep=''),
                     paste('StdError',1:gmcpSimObj$nHypothesis,sep=''),
                     paste('TestStat',1:gmcpSimObj$nHypothesis,sep=''),
                     paste('RawPvalues',1:gmcpSimObj$nHypothesis,sep=''),
                     paste('RejStatus',1:gmcpSimObj$nHypothesis,sep=''))
  SummStatBlank <- data.frame(matrix(nrow = 0, ncol = length(SummStatNames)))
  names(SummStatBlank) <- SummStatNames

  ArmWiseDataName <- c('SimID','LookID','ArmID','EpID','Completers','Mean', 'SumOfSquares')
  ArmWiseDataBlank <- data.frame(matrix(nrow = 0, ncol = length(ArmWiseDataName)))
  names(ArmWiseDataBlank) <- ArmWiseDataName

  addList <- list(
    'CurrentLook'= 1,
    'HypoPresent' = rep(TRUE,gmcpSimObj$nHypothesis),
    'ArmsPresent' = rep(TRUE, gmcpSimObj$nArms),
    'ArmsRetained' = rep(FALSE, gmcpSimObj$nArms),
    'p_raw' = NA,
    'WH_Prev' = preSimObjs$WH,
    'rej_flag_Prev'= rep(FALSE, gmcpSimObj$nHypothesis),
    'rej_flag_Curr'= rep(FALSE, gmcpSimObj$nHypothesis),
    'DropedFlag'= rep(FALSE, gmcpSimObj$nHypothesis),
    'LastLook'=gmcpSimObj$nLooks,
    'Modify'=F,
    'ModificationLook'=c(),
    'newWeights' =NA,
    'newG' = NA,
    'SummStatBlank' = SummStatBlank,
    'ArmWiseDataBlank' = ArmWiseDataBlank,
    'SummStatDF'=NA,
    'ArmDataDF'=NA,
    'ContTrial' = T
  )

  mcpObj <- append(append(gmcpSimObj, preSimObjs), addList)
  mcpObj
}

#---------------  -
#Response Generation for interim looks based on the available arms or re-allocated sample size(Implicit SSR)
#--------------- -
getInterimSSIncr <- function(lookID,PlanSSIncr,ArmsPresent,ArmsRetained,
                                Arms.alloc.ratio,ImplicitSSR)
{
  if(ncol(PlanSSIncr) != length(ArmsPresent) || ncol(PlanSSIncr) != length(Arms.alloc.ratio)) stop("'ncol(PlanSSIncr)', length(Arms.alloc.ratio) and length(ArmsPresent) are not same")
  if(all(ArmsPresent == F) || all(ArmsRetained==T)) stop('No Arms present to Continue')

  planSS <- SS.arm <-  PlanSSIncr[lookID,] #If no Implicit SSR is needed SS.arm = Planned Sample Size

  if(ImplicitSSR == 'Selection' & lookID >1) #Re-allocate samples from only de-selected arms
  {
    AditionalSS <- sum(planSS[which(ArmsRetained)])
    ss_frac <- Arms.alloc.ratio[ArmsPresent]/sum(Arms.alloc.ratio[ArmsPresent])
    SS.arm.add <- round(AditionalSS*ss_frac)
    SS.arm.add[1] <- AditionalSS - sum(SS.arm.add[-1])
    SS.arm[which(ArmsPresent)] <- SS.arm[which(ArmsPresent)]+SS.arm.add
    SS.arm[which(!ArmsPresent)] <- NA
  }else if(ImplicitSSR == 'All' & lookID >1) #Allocate all the samples planed to the available arms
  {
    AditionalSS <- sum(planSS[which(!ArmsPresent)])
    ss_frac <- Arms.alloc.ratio[ArmsPresent]/sum(Arms.alloc.ratio[ArmsPresent])
    SS.arm.add <- round(AditionalSS*ss_frac)
    SS.arm.add[1] <- AditionalSS - sum(SS.arm.add[-1])
    SS.arm[which(ArmsPresent)] <- SS.arm[which(ArmsPresent)]+SS.arm.add
    SS.arm[which(!ArmsPresent)] <- NA
  }else
  {
    SS.arm[!ArmsPresent] <- NA
  }
  SS.arm
}

ifElse <- function(test,yes,no)
{
  if(test){
    yes
  }else
  {
    no
  }
}












