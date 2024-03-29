

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
genIncrLookSummary <- function(SimSeed,
                                  simID,
                                  lookID,
                                  Arms.Mean,
                                  Arms.std.dev,
                                  Arms.Prop,
                                  Arms.alloc.ratio,
                                  Arms.SS,
                                  EPCorr,
                                  ArmsPresent,
                                  HypoPresent,
                                  HypoMap)
{
  tryCatch(
    {
      returnSubjData <- F #Flag to get the Subject level responses(Arm-Wise) as a return object
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

          #The following code is commented as it is specific to all normal endpoint
          # arm_mean <- sapply(Arms.Mean[grps], function(x)x[armIDX])
          # arm_sd <- sapply(Arms.std.dev[grps], function(x)x[armIDX])
          # if(length(grps)==1)
          # {
          #   arm_Sigma <- as.numeric(arm_sd)
          # }else
          # {
          #   arm_corr <- EPCorr[grps,grps]
          #   arm_Sigma <- EPCorr[grps,grps] * (arm_sd %*% t(arm_sd))
          # }
          #-------------------------------------------------------------------

          arm_SS <- as.numeric(Arms.SS[armIDX])
          # arm_response <- mvtNormResponse(n = arm_SS,
          #                                 mean = arm_mean,
          #                                 arm_sigma = arm_Sigma,
          #                                 seed = getRunSeed(SimSeed = SimSeed,simID = simID, lookID = lookID,armIndex = armIDX))
          vEPType <- sapply(grps, function(grpIDX){
            unique(HypoMapAvl_arm$EpType[HypoMapAvl_arm$Groups == grpIDX])
          })
          arm_response <- genNormToOther2(nArmID = armIDX,
                                          nSubject = arm_SS,
                                          vEPs = grps,
                                          vEPType = vEPType,
                                          lNormMean = Arms.Mean,
                                          lNormStdDev = Arms.std.dev,
                                          lProp = Arms.Prop,
                                          mNormCorr = EPCorr,
                                          nSeed = getRunSeed(SimSeed = SimSeed,simID = simID, lookID = lookID,armIndex = armIDX))


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

      #The following code is sensitive to the column names
      colnames(SoS_df) <- colnames(SUM_df) <- colnames(MEAN_df) <- c("Control", "Treatment")

      if(returnSubjData){
        return(list('SSIncrDF'=SS_df, 'SumOfSquareDF'=SoS_df, 'SumDF'=SUM_df, 'MeanDF'=MEAN_df, 'ArmData'=armData, 'SubjData'=SubjData))
      }else{
        return(list('SSIncrDF'=SS_df, 'SumOfSquareDF'=SoS_df, 'SumDF'=SUM_df, 'MeanDF'=MEAN_df, 'ArmData'=armData))
      }

    },
    error = function(err){
      print('Error in genIncrLookSummary')
      traceback()
    }
  )

}

#----------------- -
#Compute Test Statistics
#----------------- -

getPerLookTestStat <- function(simID,
                               lookID,
                               TestStatCont,
                               TestStatBin,
                               Arms.std.dev,
                               IncrLookSummary,
                               IncrLookSummaryPrev,
                               HypoMap,
                               Cumulative)
{
 tryCatch(
   {
     if(!Cumulative || lookID==1)
     {
       #delta incremental scale
       delta_incr <- (IncrLookSummary$MeanDF$Treatment - IncrLookSummary$MeanDF$Control)

       vEPs <- unique(HypoMap$Groups)
       vEpType <- sapply(vEPs, function(epIDX){
         unique(HypoMap$EpType[HypoMap$Groups == epIDX])
       })

       #The row position should be same as IncrLookSummary$ArmData
       ArmDataMap <- dplyr::right_join(IncrLookSummary$ArmData,
                                       data.frame('EpID'=vEPs, 'EpType'=vEpType),
                                       by = 'EpID')

       std_incr_arm <- nt <- c()
       #Compute arm-wise standard deviations
       for(armIDX in 1:nrow(IncrLookSummary$ArmData))
       {
         if(ArmDataMap$EpType[armIDX]=='Continuous'){
           #Continuous-------------------------------------------------
           if(TestStatCont=='z'){
             n <- as.numeric(IncrLookSummary$ArmData$Completers[armIDX])
             std_incr_arm <- c(std_incr_arm, Arms.std.dev[[IncrLookSummary$ArmData$EpID[armIDX]]]
                               [IncrLookSummary$ArmData$ArmID[armIDX]])
             nt <- c(nt, n)
           }else
           {
             n <- as.numeric(IncrLookSummary$ArmData$Completers[armIDX])
             SoS <- as.numeric(IncrLookSummary$ArmData$SumOfSquares[armIDX])
             x_bar <- as.numeric(IncrLookSummary$ArmData$Mean[armIDX])
             std_incr_arm <- c(std_incr_arm, sqrt((1/(n-1))*(SoS-n*x_bar^2)))
             nt <- c(nt, n)
           }

         }else if(ArmDataMap$EpType[armIDX]=='Binary'){
           #Binary----------------------------------------------------
           n <- as.numeric(ArmDataMap$Completers[armIDX])
           prop <- as.numeric(ArmDataMap$Mean[armIDX])
           std_incr_arm <- c(std_incr_arm, prop*(1-prop)/n)
           nt <- c(nt, n)
         }
       }

       #Return output of the following loop hIDX in 1:nrow(HypoMap)
       SE_Pair_incr <- testStatIncr <- pValIncr <- rep(NA, length(delta_incr))

       for(hIDX in 1:nrow(HypoMap)){
         if(!is.na(delta_incr[hIDX])){
           epIDX <- HypoMap$Groups[hIDX]
           epTypeIDX <- HypoMap$EpType[hIDX]
           ctrIDX <- HypoMap$Control[hIDX]
           trtIDX <- HypoMap$Treatment[hIDX]

           if(epTypeIDX == 'Continuous')
           {
             #Test Stat for Continuous Endpoints----------------------------
             if(TestStatCont == 'z' || TestStatCont == 't-unequal'){
               sigma_ctr <- std_incr_arm[which(IncrLookSummary$ArmData$EpID == epIDX &
                                                 IncrLookSummary$ArmData$ArmID == ctrIDX)]
               n_ctr <- IncrLookSummary$ArmData$Completers[which(IncrLookSummary$ArmData$EpID == epIDX &
                                                                   IncrLookSummary$ArmData$ArmID == ctrIDX)]

               n_trt <- IncrLookSummary$ArmData$Completers[which(IncrLookSummary$ArmData$EpID == epIDX &
                                                                   IncrLookSummary$ArmData$ArmID == trtIDX)]
               sigma_trt <- std_incr_arm[which(IncrLookSummary$ArmData$EpID == epIDX &
                                                 IncrLookSummary$ArmData$ArmID == trtIDX)]

               SE <- sqrt(((sigma_trt^2)/n_trt + (sigma_ctr^2)/n_ctr))
               testStat <- delta_incr[hIDX]/SE

               SE_Pair_incr[hIDX] <- SE
               testStatIncr[hIDX] <-  testStat

               if(TestStatCont == 'z') #Z-stat
               {
                 pVal <- 1-pnorm(testStat)

               }else #t unequal
               {
                 V <- ((sigma_trt^2)/n_trt + (sigma_ctr^2)/n_ctr)
                 b <- (((sigma_trt^2)/n_trt)^2)/(n_trt-1) +
                   (((sigma_ctr^2)/n_ctr)^2)/(n_ctr-1)
                 dof <- round(V^2/b)
                 pVal <- 1-pt(q = testStat, df = dof)
               }
               pValIncr[hIDX] <- pVal
               #End of computations for Z and t-unequal

             }else if(TestStatCont == 't-equal')
             {
               sigma_ep <- std_incr_arm[which(IncrLookSummary$ArmData$EpID == epIDX)]
               n_ep <- nt[which(IncrLookSummary$ArmData$EpID == epIDX)]
               sigma_pool <- sqrt(sum(((n_ep-1)*sigma_ep^2))/sum(n_ep-1))

               n_ctr <- IncrLookSummary$ArmData$Completers[which(IncrLookSummary$ArmData$EpID == epIDX &
                                                                   IncrLookSummary$ArmData$ArmID == ctrIDX)]

               n_trt <- IncrLookSummary$ArmData$Completers[which(IncrLookSummary$ArmData$EpID == epIDX &
                                                                   IncrLookSummary$ArmData$ArmID == trtIDX)]

               SE <- sqrt(1/n_trt + 1/n_ctr)*sigma_pool
               testStat <- delta_incr[hIDX]/SE

               SE_Pair_incr[hIDX] <- SE
               testStatIncr[hIDX] <-  testStat

               dof <- sum(n_ep)-length(nt) #(n-k-1)
               pVal <- 1-pt(q = testStat, df = dof)
               pValIncr[hIDX] <- pVal
               #End of computations for t-equal

             }else
             {
               stop("Error: TestStat computation(Incr.)")
             }
             #End of epTypeIDX == 'Continuous'

           }else if(epTypeIDX == 'Binary')
           {
             if(TestStatBin == 'UnPooled'){
               #Eqation 5a P2P3 technical documentation
               sigma_ctr <- std_incr_arm[which(IncrLookSummary$ArmData$EpID == epIDX &
                                                 IncrLookSummary$ArmData$ArmID == ctrIDX)]
               sigma_trt <- std_incr_arm[which(IncrLookSummary$ArmData$EpID == epIDX &
                                                 IncrLookSummary$ArmData$ArmID == trtIDX)]
               SE <- sqrt(sigma_ctr + sigma_trt)
               testStat <- delta_incr[hIDX]/SE
               pVal <- 1-pnorm(testStat)

               SE_Pair_incr[hIDX] <- SE
               testStatIncr[hIDX] <-  testStat
               pValIncr[hIDX] <- pVal

             }else if(TestStatBin == 'Pooled'){
               #Eqation 3a P2P3 technical documentation
               prop_ctr <- IncrLookSummary$ArmData$Mean[IncrLookSummary$ArmData$EpID == epIDX &
                                                          IncrLookSummary$ArmData$ArmID == ctrIDX]
               n_ctr <- IncrLookSummary$ArmData$Completers[IncrLookSummary$ArmData$EpID == epIDX &
                                                       IncrLookSummary$ArmData$ArmID == ctrIDX]

               prop_trt <- IncrLookSummary$ArmData$Mean[IncrLookSummary$ArmData$EpID == epIDX &
                                                          IncrLookSummary$ArmData$ArmID == trtIDX]
               n_trt <- IncrLookSummary$ArmData$Completers[IncrLookSummary$ArmData$EpID == epIDX &
                                                       IncrLookSummary$ArmData$ArmID == trtIDX]

               pooledProp <- (prop_ctr*n_ctr + prop_trt*n_trt)/(n_ctr + n_trt)

               SE <- sqrt(pooledProp*(1-pooledProp)*(1/n_ctr + 1/n_trt))

               testStat <- delta_incr[hIDX]/SE
               pVal <- 1-pnorm(testStat)

               SE_Pair_incr[hIDX] <- SE
               testStatIncr[hIDX] <-  testStat
               pValIncr[hIDX] <- pVal
             }
             #End of epTypeIDX == 'Binary'
           }
         }
       }

       # std_err_incr <- lapply(1:nrow(IncrLookSummary$MeanDF), function(h){
       #   #print(h)
       #   unlist(lapply(1:2, function(trtIDX){
       #     #print(trtIDX)
       #     sqrt((1/(as.numeric(IncrLookSummary$SSIncrDF[h,trtIDX])-1))*
       #            (as.numeric(IncrLookSummary$SumOfSquareDF[h,trtIDX])-
       #               as.numeric(IncrLookSummary$SSIncrDF[h,trtIDX])*
       #               (as.numeric(IncrLookSummary$MeanDF[h,trtIDX])^2)))
       #   }))})
       #
       # #Pair-Wise Independent Std.dev for test stat
       # SE_Pair_incr <- unlist(lapply(1:nrow(IncrLookSummary$MeanDF), function(h){
       #   sqrt(std_err_incr[[h]][1]^2/as.numeric(IncrLookSummary$SSIncrDF[h,1])+
       #          std_err_incr[[h]][2]^2/as.numeric(IncrLookSummary$SSIncrDF[h,2]))}))
       #
       # #Test Stat
       # testStatIncr <- delta_incr/SE_Pair_incr
       #
       # #P-value
       # RightTail=T
       # if(RightTail)
       # {
       #   pValIncr <- 1-pnorm(testStatIncr)
       # }

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
       ctr_ss_cum <- unlist(IncrLookSummary$SSIncrDF$CtrSS)+
         unlist(IncrLookSummaryPrev$SSIncrDF$CtrSS)
       trt_ss_cum <- unlist(IncrLookSummary$SSIncrDF$TrtSS)+
         unlist(IncrLookSummaryPrev$SSIncrDF$TrtSS)
       cum_ss <- data.frame(cbind(ctr_ss_cum,trt_ss_cum))

       #Cumulative mean
       ctr_mean_cum <-  (IncrLookSummary$SumDF$Control+IncrLookSummaryPrev$SumDF$Control)/ctr_ss_cum
       trt_mean_cum <-  (IncrLookSummary$SumDF$Treatment+IncrLookSummaryPrev$SumDF$Treatment)/trt_ss_cum
       cum_mean <- data.frame(cbind(ctr_mean_cum,trt_mean_cum))

       #Cumulative delta
       delta_cum <- (trt_mean_cum-ctr_mean_cum)

       #Cumulative Arm Data
       ArmDataCum <- IncrLookSummary$ArmData
       for (i in 1:nrow(ArmDataCum)) {
         armIDX <-  ArmDataCum[i,'ArmID']; epIDX <- ArmDataCum[i,'EpID']
         rowIDX <-  which(IncrLookSummaryPrev$ArmData$ArmID == armIDX &
                      IncrLookSummaryPrev$ArmData$EpID == epIDX)
         ArmDataCum[i,'Completers'] <- IncrLookSummary$ArmData[i,'Completers']+
                                        IncrLookSummaryPrev$ArmData[rowIDX,'Completers']

         ArmDataCum[i,'Mean'] <- (IncrLookSummary$ArmData[i,'Completers']*IncrLookSummary$ArmData[i,'Mean']+
                                  IncrLookSummaryPrev$ArmData[rowIDX,'Completers']*IncrLookSummaryPrev$ArmData[rowIDX,'Mean'])/ArmDataCum[i,'Completers']

         ArmDataCum[i,'SumOfSquares'] <- IncrLookSummary$ArmData[i,'SumOfSquares']+
           IncrLookSummaryPrev$ArmData[rowIDX,'SumOfSquares']

       }

       vEPs <- unique(HypoMap$Groups)
       vEpType <- sapply(vEPs, function(epIDX){
         unique(HypoMap$EpType[HypoMap$Groups == epIDX])
       })

       #The row position should be same as IncrLookSummary$ArmData
       ArmDataMap <- dplyr::right_join(ArmDataCum,
                                       data.frame('EpID'=vEPs, 'EpType'=vEpType),
                                       by = 'EpID')

       std_cum_arm <- nt <- c()
       for(armIDX in 1:nrow(ArmDataCum))
       {
         if(ArmDataMap$EpType[armIDX]=='Continuous'){
           #Continuous-------------------------------------------------
           if(TestStatCont=='z'){
             n <- as.numeric(ArmDataCum$Completers[armIDX])
             std_cum_arm <- c(std_cum_arm, Arms.std.dev[[ArmDataCum$EpID[armIDX]]]
                              [ArmDataCum$ArmID[armIDX]])
             nt <- c(nt, n)
           }else
           {
             n <- as.numeric(ArmDataCum$Completers[armIDX])
             SoS <- as.numeric(ArmDataCum$SumOfSquares[armIDX])
             x_bar <- as.numeric(ArmDataCum$Mean[armIDX])
             std_cum_arm <- c(std_cum_arm, sqrt((1/(n-1))*(SoS-n*x_bar^2)))
             nt <- c(nt, n)
           }
         }else if(ArmDataMap$EpType[armIDX]=='Binary'){
           #Binary----------------------------------------------------
           n <- as.numeric(ArmDataMap$Completers[armIDX])
           prop <- as.numeric(ArmDataMap$Mean[armIDX])
           std_cum_arm <- c(std_cum_arm, prop*(1-prop)/n)
           nt <- c(nt, n)
         }

       }

       #Return output of the following loop hIDX in 1:nrow(HypoMap)
       SE_Pair_cum <- test_Stat_cum <- pVal_cum <- rep(NA, length(delta_cum))

       for(hIDX in 1:nrow(HypoMap)){
         if(!is.na(delta_cum[hIDX])){
           epIDX <- HypoMap$Groups[hIDX]
           epTypeIDX <- HypoMap$EpType[hIDX]
           ctrIDX <- HypoMap$Control[hIDX]
           trtIDX <- HypoMap$Treatment[hIDX]

           if(epTypeIDX == 'Continuous'){
             #Test Stat for Continuous Endpoints----------------------------
             if(TestStatCont == 'z' || TestStatCont == 't-unequal'){
               sigma_ctr <- std_cum_arm[which(ArmDataCum$EpID == epIDX &
                                                ArmDataCum$ArmID == ctrIDX)]
               n_ctr <- ArmDataCum$Completers[which(ArmDataCum$EpID == epIDX &
                                                      ArmDataCum$ArmID == ctrIDX)]

               n_trt <- ArmDataCum$Completers[which(ArmDataCum$EpID == epIDX &
                                                      ArmDataCum$ArmID == trtIDX)]
               sigma_trt <- std_cum_arm[which(ArmDataCum$EpID == epIDX &
                                                ArmDataCum$ArmID == trtIDX)]

               SE <- sqrt(((sigma_trt^2)/n_trt + (sigma_ctr^2)/n_ctr))
               testStat <- delta_cum[hIDX]/SE

               SE_Pair_cum[hIDX] <- SE
               test_Stat_cum[hIDX] <-  testStat

               if(TestStatCont == 'z')
               {
                 pVal <- 1-pnorm(testStat)

               }else #t unequal
               {
                 V <- ((sigma_trt^2)/n_trt + (sigma_ctr^2)/n_ctr)
                 b <- (((sigma_trt^2)/n_trt)^2)/(n_trt-1) +
                   (((sigma_ctr^2)/n_ctr)^2)/(n_ctr-1)
                 dof <- round(V^2/b)
                 pVal <- 1-pt(q = testStat, df = dof)
               }
               pVal_cum[hIDX] <- pVal
               #End of computations for Z and t-unequal

             }else if(TestStatCont == 't-equal')
             {
               sigma_ep <- std_cum_arm[which(ArmDataCum$EpID == epIDX)]
               n_ep <- nt[which(ArmDataCum$EpID == epIDX)]
               sigma_pool <- sqrt(sum(((n_ep-1)*sigma_ep^2))/sum(n_ep-1))

               n_ctr <- ArmDataCum$Completers[which(ArmDataCum$EpID == epIDX &
                                                      ArmDataCum$ArmID == ctrIDX)]

               n_trt <- ArmDataCum$Completers[which(ArmDataCum$EpID == epIDX &
                                                      ArmDataCum$ArmID == trtIDX)]

               SE <- sqrt(1/n_trt + 1/n_ctr)*sigma_pool
               testStat <- delta_cum[hIDX]/SE

               SE_Pair_cum[hIDX] <- SE
               test_Stat_cum[hIDX] <-  testStat

               dof <- sum(n_ep-1)
               pVal <- 1-pt(q = testStat, df = dof)
               pVal_cum[hIDX] <- pVal
               #End of computations for t-equal

             }else
             {
               stop("Error: TestStat computation(Cum.)")
               traceback()
             }
             #End of epTypeIDX == 'Continuous'
           }else if(epTypeIDX == 'Binary'){
             #Test Stat for Binary Endpoints----------------------------

             if(TestStatBin == 'UnPooled'){
               #Eqation 5a P2P3 technical documentation
               sigma_ctr <- std_cum_arm[which(ArmDataCum$EpID == epIDX &
                                                 ArmDataCum$ArmID == ctrIDX)]
               sigma_trt <- std_cum_arm[which(ArmDataCum$EpID == epIDX &
                                                 ArmDataCum$ArmID == trtIDX)]
               SE <- sqrt(sigma_ctr + sigma_trt)
               testStat <- delta_cum[hIDX]/SE
               pVal <- 1-pnorm(testStat)

               SE_Pair_cum[hIDX] <- SE
               test_Stat_cum[hIDX] <-  testStat
               pVal_cum[hIDX] <- pVal

             }else if(TestStatBin == 'Pooled'){
               #Eqation 3a P2P3 technical documentation
               prop_ctr <- ArmDataCum$Mean[ArmDataCum$EpID == epIDX &
                                                          IncrLookSummary$ArmData$ArmID == ctrIDX]
               n_ctr <- ArmDataCum$Completers[ArmDataCum$EpID == epIDX &
                                                ArmDataCum$ArmID == ctrIDX]

               prop_trt <- ArmDataCum$Mean[ArmDataCum$EpID == epIDX &
                                             ArmDataCum$ArmID == trtIDX]
               n_trt <- ArmDataCum$Completers[ArmDataCum$EpID == epIDX &
                                                ArmDataCum$ArmID == trtIDX]

               pooledProp <- (prop_ctr*n_ctr + prop_trt*n_trt)/(n_ctr + n_trt)

               SE <- sqrt(pooledProp*(1-pooledProp)*(1/n_ctr + 1/n_trt))

               testStat <- delta_cum[hIDX]/SE
               pVal <- 1-pnorm(testStat)

               SE_Pair_cum[hIDX] <- SE
               test_Stat_cum[hIDX] <-  testStat
               pVal_cum[hIDX] <- pVal
             }
             #End of epTypeIDX == 'Binary'
           }
         }
       }


       # #Cumulative sum of squares
       # ctr_SoS_cum <-  (IncrLookSummary$SumOfSquareDF$Control+IncrLookSummaryPrev$SumOfSquareDF$Control)
       # trt_SoS_cum <-  (IncrLookSummary$SumOfSquareDF$Treatment+IncrLookSummaryPrev$SumOfSquareDF$Treatment)
       # SoS_cum <- data.frame(cbind(ctr_SoS_cum,trt_SoS_cum))
       #
       # #Per arm Std.dev cumulative
       # std_err_cum <- lapply(1:nrow(cum_ss), function(h){
       #   unlist(lapply(1:2, function(trtIDX){
       #     sqrt((1/(cum_ss[h,trtIDX]-1))*
       #            (SoS_cum[h,trtIDX]-cum_ss[h,trtIDX]*(cum_mean[h,trtIDX]^2)))
       #   }))})
       #
       # #Pair-Wise Independent Std.dev(Cum.) for test stat
       # SE_Pair_cum <- unlist(lapply(1:nrow(cum_ss), function(h){
       #   sqrt(std_err_cum[[h]][1]^2/cum_ss[h,1]+std_err_cum[[h]][2]^2/cum_ss[h,2])}))
       #
       # #Test Stat cumulative
       # testStatCum <- delta_cum/SE_Pair_cum
       #
       # #P-value cumulative
       # RightTail=T
       # if(RightTail)
       # {
       #   pValCum <- 1-pnorm(testStatCum)
       # }

       sumstatdf <- data.frame(matrix(c(simID,lookID,delta_cum,SE_Pair_cum,
                                        test_Stat_cum, pVal_cum),
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
     traceback()
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

