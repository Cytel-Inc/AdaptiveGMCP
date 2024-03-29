
# The file contains supporting functions for Adaptive GMCP Simulation ----
## Author: Ajoy.M

#------------------ -
#Perform iterative simulation for Multi-Arm designs
#------------------ -
MAMSMEP_sim2 <- function(gmcpSimObj)
{
  # library(foreach)
  # library(doParallel)
  ########################Computation of intermediate inputs#####################
  preSimObjs <<- getPreSimObjs(gmcpSimObj = gmcpSimObj)

  #-----------------------------------------------------------------------------

  ########################### Run Simulations ######################################
  starttime <- Sys.time()
  SummaryStatFile <- ArmWiseSummary<- PowerTab <- data.frame()
  powersName <- c('simID','nG','nC','nD','nF')
  PowerTab <- data.frame(matrix(nrow = 0, ncol = length(powersName)))
  EfficacyTable <- data.frame()
  names(PowerTab) <- powersName
  SelectionTab <- data.frame()

  if(gmcpSimObj$Method=='CER')
  {
    if(gmcpSimObj$Parallel){
      cores=parallel::detectCores()
      cl <- parallel::makeCluster(cores[1]-1, type = "PSOCK", outfile = 'log.txt') # leave 1 core and use the remaining for this computation
      #cl <- parallel::makeCluster(clust_cores, type = "PSOCK")
      #parallel::clusterSetRNGStream(cl, iseed = preSimObjs$SimSeed)
      parallel::clusterExport(cl, c("gmcpSimObj", "preSimObjs"))

      out <- parallel::parLapply(cl = cl, 1:gmcpSimObj$nSimulation, function(x){
        out_SingleSim <- SingleSimCER(x, gmcpSimObj,preSimObjs)
        return(out_SingleSim)
      })
      parallel::stopCluster(cl)
    }else{
      out<-lapply(1:gmcpSimObj$nSimulation, function(x){
        SingleSimCER(x, gmcpSimObj,preSimObjs)})
    }
  }else
  {
    if(gmcpSimObj$Parallel){
      cores=parallel::detectCores()
      cl <- parallel::makeCluster(cores[1]-1, type = "PSOCK", outfile = 'log.txt') # leave 1 core and use the remaining for this computation
      #cl <- parallel::makeCluster(clust_cores)
      #parallel::clusterSetRNGStream(cl, iseed = preSimObjs$SimSeed)
      parallel::clusterExport(cl, c("gmcpSimObj", "preSimObjs"))

      out <- parallel::parLapply(cl = cl, 1:gmcpSimObj$nSimulation,
        function(x){
        out_SingleSim <- SingleSimCombPValue(x, gmcpSimObj,preSimObjs)
        return(out_SingleSim)
      })
      parallel::stopCluster(cl)
    }else{
      out<-lapply(1:gmcpSimObj$nSimulation, function(x){
        SingleSimCombPValue(x, gmcpSimObj,preSimObjs)})
    }
  }
  #------------------------------------------------------------------------------

  ########################## Preparing detailed output Tables ####################
  for(i in 1:length(out))
  {
    SummaryStatFile <- rbind(SummaryStatFile, out[[i]]$SummStatDF)
    ArmWiseSummary <- rbind(ArmWiseSummary, out[[i]]$ArmWiseDF)
    PowerTab <- rbind(PowerTab, out[[i]]$powerCountDF)
    EfficacyTable <- plyr::rbind.fill(EfficacyTable,out[[i]]$EfficacyTable)
    SelectionTab <- rbind(SelectionTab,out[[i]]$SelectionDF)
  }

  #Overall Powers
  Sim_power <- SimPowers(nSimulation=gmcpSimObj$nSimulation,PowerTab=PowerTab)
  Sim_power_df <- Sim_power
  Sim_power <- knitr::kable(Sim_power, align = 'c')


  #Detailed Efficacy Table
  eff_count <- colSums(EfficacyTable[,-1])
  EffTab <- data.frame('Hypothesis'=names(eff_count),
                          'Count'=eff_count,
                          'Percentage'=100*(eff_count/nrow(EfficacyTable)))
  rownames(EffTab) <- NULL
  EffTab <- knitr::kable(EffTab, align = 'c')

  if(gmcpSimObj$Selection){
    #Selection Table
    SelcCount <- table(SelectionTab$SelectedHypothesis)
    SelcPerc <- 100*(SelcCount/gmcpSimObj$nSimulation)
    SelecTab <- data.frame('Hypothesis'=names(SelcCount),
                           'Count'= as.vector(SelcCount),
                           'Percentage'= as.vector(SelcPerc))

    SelecTab <- knitr::kable(SelecTab, align = 'c')
  }else
  {
    SelecTab <- NA
  }



  elapsedTime <- Sys.time()-starttime
  #------------------------------------------------------------------------------

  if(gmcpSimObj$Method=='CER')
  {
    detailOutput <- ifElse(gmcpSimObj$SummaryStat,
                           list(
                             'PlanSampleSizeCum' = preSimObjs$planSS$CumulativeSamples,
                             'PlannedSigma' = preSimObjs$Sigma,
                             'Boundary_Table' = preSimObjs$plan_Bdry$PlanBdryTab,
                             'Overall_Powers' = Sim_power,
                             'Overall_Powers_df' = Sim_power_df,
                             'EfficacyTable'=EffTab,
                             'SelectionTable'= SelecTab,
                             'Summary_Stat' = SummaryStatFile,
                             'ArmWiseSummary' = ArmWiseSummary,
                             'Seed'= preSimObjs$SimSeed,
                             'elapsedTime' = elapsedTime
                           ),
                           list(
                             'PlanSampleSizeCum' = preSimObjs$planSS$CumulativeSamples,
                             'PlannedSigma' = preSimObjs$Sigma,
                             'Boundary_Table' = preSimObjs$plan_Bdry$PlanBdryTab,
                             'Overall_Powers' = Sim_power,
                             'Overall_Powers_df' = Sim_power_df,
                             'EfficacyTable'=EffTab,
                             'SelectionTable'= SelecTab,
                             'Seed'= preSimObjs$SimSeed,
                             'elapsedTime' = elapsedTime
                           )
    )
  }else
  {
    detailOutput <- ifElse(gmcpSimObj$SummaryStat,
                           list(
                             'PlanSampleSizeIncr' = preSimObjs$planSS$IncrementalSamples,
                             'PlannedCorrelation' = preSimObjs$PlanCorrelation,
                             'Boundary_Table' = preSimObjs$pValBdry$pValueBdryTab,
                             'Inverse_Normal_Weights' = preSimObjs$InvNormWeights$InvNormWeightsTab,
                             'Overall_Powers' = Sim_power,
                             'Overall_Powers_df' = Sim_power_df,
                             'EfficacyTable'=EffTab,
                             'SelectionTable'= SelecTab,
                             'Summary_Stat' = SummaryStatFile,
                             'ArmWiseSummary' = ArmWiseSummary,
                             'Seed'= preSimObjs$SimSeed,
                             'elapsedTime' = elapsedTime
                           ),
                           list(
                             'PlanSampleSize' = preSimObjs$planSS$IncrementalSamples,
                             'PlannedCorrelation' = preSimObjs$PlanCorrelation,
                             'Boundary_Table' = preSimObjs$pValBdry$pValueBdryTab,
                             'Inverse_Normal_Weights' = preSimObjs$InvNormWeights$InvNormWeightsTab,
                             'Overall_Powers' = Sim_power,
                             'Overall_Powers_df' = Sim_power_df,
                             'EfficacyTable'=EffTab,
                             'SelectionTable'= SelecTab,
                             'Seed'= preSimObjs$SimSeed,
                             'elapsedTime' = elapsedTime
                           )
    )
  }
  if(gmcpSimObj$plotGraphs)
  {
    list('DetailOutTabs'=detailOutput, 'iniGraph' = preSimObjs$iniGraph)
  }else
  {
    list('DetailOutTabs'=detailOutput)
  }
}



getPreSimObjs <- function(gmcpSimObj)
{
  #----------------------------------------------------------------------------
  ################## Generate unique Seed for the run #########################
  SimSeed <- ifelse(gmcpSimObj$Seed=='Random',sample(1:10000,1),gmcpSimObj$Seed)

  #----------------------------------------------------------------------------
  ################### Arms-EPs-Hypotheses Mapping #####################################
  HypoMap <- getHypoMap2(des.type = gmcpSimObj$des.type,
             nHypothesis = gmcpSimObj$nHypothesis,
             nEps = gmcpSimObj$nEps,
             nArms = gmcpSimObj$nArms,
             lEpType = gmcpSimObj$lEpType)


  #------------------------------------------------------------------------------
  #################### Identify True Null based on given Response #######################
  #TrueNull <- checkTrueNull2(HypoMap = HypoMap, Arms.Mean = gmcpSimObj$Arms.Mean)
  TrueNull <- checkTrueNull3(HypoMap = HypoMap,
                             Arms.Mean = gmcpSimObj$Arms.Mean,
                             Arms.Prop = gmcpSimObj$Arms.Prop)


  #---------------------------------------------------------------------------------------
  #################### Weights for all intersection hypothesis ##########################
  ############# Weights for all intersection hypothesis#########################
  allGraphs <- genWeights(w = gmcpSimObj$IntialWeights,
                          g = gmcpSimObj$G,
                          HypothesisName = HypoMap$Hypothesis)
  WH <- allGraphs$IntersectionWeights

  #WH <- as.data.frame(gMCPLite::generateWeights(g = gmcpSimObj$G,w = gmcpSimObj$IntialWeights))
  #colnames(WH) <- c(paste("H",1:gmcpSimObj$nHypothesis, sep=''),
  #                  paste("Weight",1:gmcpSimObj$nHypothesis, sep = ''))

  #--------------------------------------------------------------------------------------
  ############################## Planned Sample Size ################################
  planSS <- getPlanAllocatedSamples(SS = gmcpSimObj$Max_SS,
                                    allocRatio = gmcpSimObj$Arms.alloc.ratio,
                                    info_frac = gmcpSimObj$InfoFrac)

  planSSIncr <- planSS$IncrementalSamples
  planSSCum <- planSS$CumulativeSamples

  #----------------------------------------------------------------------------------
  ########Computation of Inverse Normal Weights from planned sample size##############
  #Needed for combining p-values method and CER with FWER control as 'CombinationTest'
  InvNormWeights <- getInvNormWeights(planSSIncr = planSSIncr)

  #----------------------------------------------------------------------------------

  #Prepare the preSim objects for two methods
  if(gmcpSimObj$Method == 'CombPValue') #Inputs for Combining P-values method
  {
    #----------------------------------------------------------------------------------
    ######################Get the stage-wise p-value boundaries####################
    pValBdry <- getPvalBdry(alpha = gmcpSimObj$alpha,
                nLooks = gmcpSimObj$nLooks,
                info_frac = gmcpSimObj$InfoFrac,
                typeOfDesign = gmcpSimObj$typeOfDesign)

    #----------------------------------------------------------------------------------
    ######################Get the stage-wise incremental Correlation####################
    PlanCorrelation <- getPlanCorrelation(nHypothesis = gmcpSimObj$nHypothesis,
                                          SS_Incr = planSSIncr,
                                          Arms.std.dev = gmcpSimObj$Arms.std.dev,
                                          test.type = gmcpSimObj$test.type,
                                          EpType = gmcpSimObj$lEpType,
                                          prop.ctr = gmcpSimObj$prop.ctr)

    #----------------------------------------------------------------------------------
    PreSimObj <- list(
      'IndexSet' = HypoMap$Hypothesis,
      'SimSeed' = SimSeed,    'HypoMap' = HypoMap,
      'TrueNull' = TrueNull,  'planSS' = planSS,
      'WH' = WH,              'InvNormWeights' = InvNormWeights,
      'pValBdry'=pValBdry,     'PlanCorrelation' = PlanCorrelation,
      'W_Norm'=InvNormWeights$W_Norm,
      'AdjPValues'=NA
    )

  }else if(gmcpSimObj$Method == 'CER') #Inputs for CER method
  {
    #----------------------------------------------------------------------------------
    ########Computation of covariance matrix##############
    if(gmcpSimObj$test.type == 'Partly-Parametric' || gmcpSimObj$test.type == 'Parametric')
    {
      Sigma <- getSigma(SS_Cum = planSS$CumulativeSamples,
                        sigma = gmcpSimObj$Arms.std.dev,
                        allocRatio = gmcpSimObj$Arms.alloc.ratio,
                        EpType = gmcpSimObj$lEpType,
                        prop.ctr = gmcpSimObj$prop.ctr)

    }else
    {
      Sigma <- NA
    }

    #----------------------------------------------------------------------------------
    ########Computation of Planned Boundaries##############
    plan_Bdry <- planBdryCER(nHypothesis = gmcpSimObj$nHypothesis,
                             nEps=gmcpSimObj$nEps,
                             nLooks=gmcpSimObj$nLooks,
                             alpha=gmcpSimObj$alpha,
                             info_frac=gmcpSimObj$InfoFrac,
                             typeOfDesign=gmcpSimObj$typeOfDesign,
                             test.type=gmcpSimObj$test.type,
                             Sigma=Sigma,
                             WH=WH,
                             HypoMap=HypoMap,
                             Scale = 'Score')



    #----------------------------------------------------------------------------------
    PreSimObj <- list(
      'IndexSet' = HypoMap$Hypothesis,
      'IntialHypothesis'= HypoMap$Hypothesis,
      'SimSeed' = SimSeed,                       'HypoMap' = HypoMap,
      'TrueNull' = TrueNull,                     'planSS' = planSS,
      'WH' = WH,                                 'Sigma'=Sigma,
      'plan_Bdry'=plan_Bdry,                     'InvNormWeights' = InvNormWeights
    )

  }

  #Plot the initial Graph
  SubText <- getPlotText(HypoMap)
  if(gmcpSimObj$plotGraphs)
  {
    iniGraph <- plotGraph(HypothesisName = HypoMap$Hypothesis,
              w = gmcpSimObj$IntialWeights,
              G = gmcpSimObj$G,
              Titel = 'Initial Graph',
              Text = SubText)
    PreSimObj$iniGraph <- iniGraph
  }

  PreSimObj
}











