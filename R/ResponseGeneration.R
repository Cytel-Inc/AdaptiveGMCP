

#----------------- -
#Generate Arm-Wise Multivariate Random Normal responses
#----------------- -
mvtNormResponse <- function(n, mean, arm_sigma, seed)
{
  #Make Seed restricted to local##
  # old <- .Random.seed
  # on.exit( { .Random.seed <<- old } )
  set.seed(seed = seed)
  #################################
  mvtnorm::rmvnorm(n = n, mean = as.numeric(mean), sigma = as.matrix(arm_sigma))
}

#----------------- -
# Arm-Wise summary
#----------------- -
armSumry <- function(x)
{
  arm_sum <- apply(x, 2, sum)
  arm_avg <- apply(x, 2, mean)
  arm_sumOfsquare <-  apply(x^2, 2, sum)
  list('Avg'=arm_avg, 'Sum'=arm_sum, 'SumOfSquare'=arm_sumOfsquare)
}


#----------------- -
# Arm-Wise summary ArmData
#----------------- -
genIncrLookSummaryDOM <- function(SimSeed,simID,lookID,Arms.Mean,Arms.std.dev,Arms.alloc.ratio,Arms.SS,
                              EPCorr,ArmsPresent,HypoPresent,HypoMap)
{
  tryCatch(
    {
      returnSubjData <- T
      #filter by available hypothesis
      HypoMapAvl <- HypoMap[HypoPresent,]

      CtrSS <- TrtSS <- CtrSoS <- TrtSoS <- CtrSum <- TrtSum <-
        CtrMean <- TrtMean <- rep(NA, nrow(HypoMap))

      armData <- data.frame()
      SubjData <- list()
      for(armIDX in 1:length(ArmsPresent))
      {
        if(ArmsPresent[armIDX])
        {
          #Filter by arm index
          HypoMapAvl_arm <- ifElse(armIDX %in% unique(HypoMap$Control),
                                   HypoMapAvl[HypoMapAvl$Control==armIDX,],
                                   HypoMapAvl[HypoMapAvl$Treatment==armIDX,])

          grps <- unique(HypoMapAvl_arm$Groups)

          arm_mean <- sapply(Arms.Mean[grps], function(x)x[armIDX])
          arm_sd <- sapply(Arms.std.dev[grps], function(x)x[armIDX])
          if(length(grps)==1)
          {
            arm_Sigma <- as.numeric(arm_sd)
          }else
          {
            arm_corr <- EPCorr[grps,grps]
            arm_Sigma <- EPCorr[grps,grps] * (arm_sd %*% t(arm_sd))
          }

          arm_SS <- as.numeric(Arms.SS[armIDX])
          arm_response <- mvtNormResponse(n = arm_SS,
                                          mean = arm_mean,
                                          arm_sigma = arm_Sigma,
                                          seed = getRunSeed(SimSeed = SimSeed,simID = simID, lookID = lookID,armIndex = armIDX))
          if(returnSubjData) SubjData[[paste('Arm',armIDX,sep = '')]] <- arm_response

          arm_sumry <- data.frame(armSumry(arm_response))
          arm_sumry$Groups <- grps

          armdatadf <- data.frame(cbind(rep(simID,length(grps)), rep(lookID,length(grps)),
                                        rep(armIDX,length(grps)),grps,
                                        rep(arm_SS, length(grps)), arm_sumry$Avg,
                                        arm_sumry$SumOfSquare))
          colnames(armdatadf) <- c('SimID','LookID','ArmID','EpID','Completers','Mean','SumOfSquares')
          armData <- rbind(armData,armdatadf)

          if(armIDX==1)
          {
            CtrSS[which(HypoMap$Control == armIDX)] = Arms.SS[armIDX]
          }else
          {
            TrtSS[which(HypoMap$Treatment == armIDX)] = Arms.SS[armIDX]
          }

          for (g in grps) {
            if(armIDX==1)
            {
              CtrSoS[which(HypoMap$Groups == g & HypoMap$Control == armIDX)]= arm_sumry$SumOfSquare[arm_sumry$Groups==g]
              CtrSum[which(HypoMap$Groups == g & HypoMap$Control == armIDX)]= arm_sumry$Sum[arm_sumry$Groups==g]
              CtrMean[which(HypoMap$Groups == g & HypoMap$Control == armIDX)]= arm_sumry$Avg[arm_sumry$Groups==g]

            }else
            {
              TrtSoS[which(HypoMap$Groups == g & HypoMap$Treatment == armIDX)]= arm_sumry$SumOfSquare[arm_sumry$Groups==g]
              TrtSum[which(HypoMap$Groups == g & HypoMap$Treatment == armIDX)]= arm_sumry$Sum[arm_sumry$Groups==g]
              TrtMean[which(HypoMap$Groups == g & HypoMap$Treatment == armIDX)]= arm_sumry$Avg[arm_sumry$Groups==g]
            }
          }
        }
        #End of Arm-wise for loop#
      }

      SS_df <- data.frame(cbind(CtrSS,TrtSS))
      SoS_df <- data.frame(cbind(CtrSoS,TrtSoS))
      SUM_df <- data.frame(cbind(CtrSum,TrtSum))
      MEAN_df <- data.frame(cbind(CtrMean,TrtMean))
      colnames(SoS_df) <- colnames(SUM_df) <- colnames(MEAN_df) <- names(HypoMap)[3:4]

      if(returnSubjData){
        return(list('SSIncrDF'=SS_df, 'SumOfSquareDF'=SoS_df, 'SumDF'=SUM_df, 'MeanDF'=MEAN_df, 'ArmData'=armData, 'SubjData'=SubjData))
      }else{
        return(list('SSIncrDF'=SS_df, 'SumOfSquareDF'=SoS_df, 'SumDF'=SUM_df, 'MeanDF'=MEAN_df, 'ArmData'=armData))
      }

    },
    error = function(err){
      print('Error in genIncrLookSummaryDOM')
    }
  )

}

#----------------- -
#Compute Test Statistics
#----------------- -

getPerLookTestStatDOM <- function(simID, lookID, IncrLookSummaryDOM,
                                  IncrLookSummaryDOMPrev,Cumulative)
{
 tryCatch(
   {
     if(!Cumulative || lookID==1)
     {
       #delta incremental scale
       delta_incr <- (IncrLookSummaryDOM$MeanDF$Treatment - IncrLookSummaryDOM$MeanDF$Control)

       #Per arm Std.dev
       std_err_incr <- lapply(1:nrow(IncrLookSummaryDOM$MeanDF), function(h){
         #print(h)
         unlist(lapply(1:2, function(trtIDX){
           #print(trtIDX)
           sqrt((1/(as.numeric(IncrLookSummaryDOM$SSIncrDF[h,trtIDX])-1))*
                  (as.numeric(IncrLookSummaryDOM$SumOfSquareDF[h,trtIDX])-
                     as.numeric(IncrLookSummaryDOM$SSIncrDF[h,trtIDX])*(as.numeric(IncrLookSummaryDOM$MeanDF[h,trtIDX])^2)))
         }))})

       #Pair-Wise Independent Std.dev for test stat
       SE_Pair_incr <- unlist(lapply(1:nrow(IncrLookSummaryDOM$MeanDF), function(h){
         sqrt(std_err_incr[[h]][1]^2/as.numeric(IncrLookSummaryDOM$SSIncrDF[h,1])+
                std_err_incr[[h]][2]^2/as.numeric(IncrLookSummaryDOM$SSIncrDF[h,2]))}))

       #Test Stat
       testStatIncr <- delta_incr/SE_Pair_incr

       #P-value
       RightTail=T
       if(RightTail)
       {
         pValIncr <- 1-pnorm(testStatIncr)
       }

       sumstatdf <- data.frame(matrix(c(simID,lookID,delta_incr,SE_Pair_incr,testStatIncr, pValIncr),
                                      nrow = 1))
       colnames(sumstatdf) <- c('SimID','LookID',
                                paste('Delta',1:length(delta_incr),sep=''),
                                paste('StdError',1:length(delta_incr),sep=''),
                                paste('TestStat',1:length(delta_incr),sep=''),
                                paste('RawPvalues',1:length(delta_incr),sep=''))


     }else #Cumulative Scale test stat
     {
       #Cumulative Sample Size
       ctr_ss_cum <- unlist(IncrLookSummaryDOM$SSIncrDF$CtrSS)+
         unlist(IncrLookSummaryDOMPrev$SSIncrDF$CtrSS)
       trt_ss_cum <- unlist(IncrLookSummaryDOM$SSIncrDF$TrtSS)+
         unlist(IncrLookSummaryDOMPrev$SSIncrDF$TrtSS)
       cum_ss <- data.frame(cbind(ctr_ss_cum,trt_ss_cum))

       #Cumulative mean
       ctr_mean_cum <-  (IncrLookSummaryDOM$SumDF$Control+IncrLookSummaryDOMPrev$SumDF$Control)/ctr_ss_cum
       trt_mean_cum <-  (IncrLookSummaryDOM$SumDF$Treatment+IncrLookSummaryDOMPrev$SumDF$Treatment)/trt_ss_cum
       cum_mean <- data.frame(cbind(ctr_mean_cum,trt_mean_cum))

       #Cumulative delta
       delta_cum <- (trt_mean_cum-ctr_mean_cum)

       #Cumulative sum of squares
       ctr_SoS_cum <-  (IncrLookSummaryDOM$SumOfSquareDF$Control+IncrLookSummaryDOMPrev$SumOfSquareDF$Control)
       trt_SoS_cum <-  (IncrLookSummaryDOM$SumOfSquareDF$Treatment+IncrLookSummaryDOMPrev$SumOfSquareDF$Treatment)
       SoS_cum <- data.frame(cbind(ctr_SoS_cum,trt_SoS_cum))

       #Per arm Std.dev cumulative
       std_err_cum <- lapply(1:nrow(cum_ss), function(h){
         unlist(lapply(1:2, function(trtIDX){
           sqrt((1/(cum_ss[h,trtIDX]-1))*
                  (SoS_cum[h,trtIDX]-cum_ss[h,trtIDX]*(cum_mean[h,trtIDX]^2)))
         }))})

       #Pair-Wise Independent Std.dev(Cum.) for test stat
       SE_Pair_cum <- unlist(lapply(1:nrow(cum_ss), function(h){
         sqrt(std_err_cum[[h]][1]^2/cum_ss[h,1]+std_err_cum[[h]][2]^2/cum_ss[h,2])}))

       #Test Stat cumulative
       testStatCum <- delta_cum/SE_Pair_cum

       #P-value cumulative
       RightTail=T
       if(RightTail)
       {
         pValCum <- 1-pnorm(testStatCum)
       }

       sumstatdf <- data.frame(matrix(c(simID,lookID,delta_cum,SE_Pair_cum,testStatCum, pValCum),
                                      nrow = 1))
       colnames(sumstatdf) <- c('SimID','LookID',
                                paste('Delta',1:length(delta_cum),sep=''),
                                paste('StdError',1:length(delta_cum),sep=''),
                                paste('TestStat',1:length(delta_cum),sep=''),
                                paste('RawPvalues',1:length(delta_cum),sep=''))


     }
     sumstatdf
   },
   error = function(err){
     print('Error in getPerLookTestStatDOM')
   }
 )
}


#----------------- -
#Summary Arm-Wise
#----------------- -
perArmData <- function(simSeed,simID,lookID,Arms.Mean,Arms.std.dev,Arms.alloc.ratio,Arms.SS,
                       ArmsPresent,HypoPresent,HypoMap)
{
  if(lookID==1) #for look 1 and FSD
  {
    SS.arm <- mcpObj$PlannedSS[1,]
    mu.arm <- gmcpSimObj$Arms.Mean
    sd.arm <- gmcpSimObj$Arms.std.dev
    result_list <- lapply(1:length(SS.arm), function(x) {normalResponse(armID = names(SS.arm)[x], n = SS.arm[x],
                                                                        mu = mu.arm[x], s = sd.arm[x],Arm.seed = getRunSeed(SimSeed, simID, lookID, x))})
    df <- data.frame(Arm = character(), Mean = numeric(), SE = numeric())
    df <- do.call(rbind, lapply(result_list, as.data.frame))
    rownames(df) <- NULL
    return(df)

  }else #interim and Final look
  {
    LookNrespose <- getInterimResposes(lookID,mcpObj,gmcpSimObj)
    mu.arm <- LookNrespose$mu.arm
    sd.arm <- LookNrespose$sd.arm
    SS.arm <- LookNrespose$SS.arm
    SS.arm <- SS.arm[names(mu.arm)]

    result_list <- lapply(1:length(SS.arm), function(x) {normalResponse(armID = names(SS.arm)[x], n = SS.arm[x],
                                                                        mu = mu.arm[x], s = sd.arm[x],Arm.seed = getRunSeed(SimSeed, simID, lookID, x))})
    df <- data.frame(Arm = character(), Mean = numeric(), SE = numeric())
    df <- do.call(rbind, lapply(result_list, as.data.frame))
    rownames(df) <- NULL
    return(df)

  }
}

