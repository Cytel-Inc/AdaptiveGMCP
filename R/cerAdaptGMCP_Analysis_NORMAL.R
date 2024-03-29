

#' Function to perform Adaptive GMCP Analysis following Combined P-Value(Inverse Normal) method
#' @param nArms Number of Arms
#' @param nEps Number of End points
#' @param SampleSize Plan Sample Size
#' @param sigma Arm-Wise sigma for each endpoint(for EpType = "Continuous")
#' @param allocRatio Arm-Wise allocation ratio
#' @param WI Vector of Initial Weights for Global Null(default = \code{rep(1/4,4)})
#' @param G  Transition Matrix (default = \code{matrix(c(0,1/3,1/3,1/3,  1/3,0,1/3,1/3, 1/3,1/3,0,1/3, 1/3,1/3,1/3,0),nrow = 4)})
#' @param test.type Character to specify the type of test want to perform; "Parametric": Weighted Dunnett , "Non-Parametric": Weighted Bonferroni and  'Partly-Parametric': Mixed type Tests.
#' @param alpha Type-1 error
#' @param info_frac Vector of information fraction
#' @param typeOfDesign The type of design. Type of design is one of the following: O'Brien & Fleming ("OF"), Pocock ("P"), Wang & Tsiatis Delta class ("WT"), Pampallona & Tsiatis ("PT"), Haybittle & Peto ("HP"), Optimum design within Wang & Tsiatis class ("WToptimum"), O'Brien & Fleming type alpha spending ("asOF"), Pocock type alpha spending ("asP"), Kim & DeMets alpha spending ("asKD"), Hwang, Shi & DeCani alpha spending ("asHSD"), no early efficacy stop ("noEarlyEfficacy"), default is "OF".
#' @param AdaptStage2 TRUE: Adaptation option will be given for stage-2, FALSE : proceed as planned.
#' @param plotGraphs TRUE: plot intermediate graphs
#' @example ./internalData/AdaptGMCP_CER_Analysis_Example.R
#' @export
adaptGMCP_CONT_CER <- function(
    nArms = 3,
    nEps  = 2,
    SampleSize = 500,
    sigma = list('EP1' = c(1,1,1), 'EP2' = c(1,1,1)),
    allocRatio = c(1,1,1),
    WI =  c(0.5,0.5,0,0),
    G = matrix(c(0,0.5,0.5,0,
                 0.5,0,0,0.5,
                 0,1,0,0,
                 1,0,0,0),
               nrow = nEps*(nArms-1), byrow = T),
    test.type = "Partly-Parametric",
    alpha = 0.025,
    info_frac = c(0.5,1),
    typeOfDesign = "asOF",
    AdaptStage2 = TRUE,
    plotGraphs = TRUE
)
{
  ###### Input Validation #####
  # stopifnot('Number of Arms must be > 2',length(nArms) <= 2)
  # stopifnot('Number of End points must be >= 1',length(nEps) < 1)

  #############################
  EpType <- "Continuous"
  TailType <- 'RightTail'       ##Default Right
  Hypothesis <- 'CommonControl' ##Default CommonControl

  #Adaptation Choice
  MultipleWinners <- TRUE
  Selection <- TRUE
  UpdateStrategy <- TRUE
  ModifySamples <- TRUE

  nLooks <- length(info_frac)
  des.type <- 'MAMSMEP'

  nHypothesis <- nEps*(nArms-1)
  GlobalIndexSet <- paste("H",1:nHypothesis, sep ='')
  ArmsPresent <- 1:nArms

  #Map for arms and hypothesis
  HypoMap <- getHypoMap(des.type=des.type, nHypothesis=nHypothesis,
                        nEps = nEps, nArms=nArms)

  #Computation of Intersection weights
  allGraphs <- genWeights(w = WI, g = G, HypothesisName=GlobalIndexSet)
  WH <- allGraphs$IntersectionWeights

  #Plot the initial Graph
  SubText <- getPlotText(HypoMap)
  if(plotGraphs)
  {
    plotGraph(HypothesisName = GlobalIndexSet,
              w = WI,
              G = G,
              Titel = 'Initial Graph',
              Text = SubText)
  }

  rej_flag_Prev <- rej_flag_Curr <- DropedFlag <- rep(FALSE, nHypothesis)
  names(rej_flag_Prev) <- names(rej_flag_Curr) <- names(DropedFlag) <-paste("H",1:nHypothesis, sep = '')

  #info to run per Stage-Wise analysis
  mcpObj <- list(
    'CurrentLook'= NA,
    'EpType' = EpType,
    'IntialHypothesis'=GlobalIndexSet,
    'test.type' = test.type,
    'IndexSet'= GlobalIndexSet,
    'ArmsPresent' = ArmsPresent,
    'sigma' = sigma,
    'allocRatio' = allocRatio,
    'Stage2allocRatio' = allocRatio,
    'p_raw' = NA,
    'WH_Prev' = WH,
    'WH'= WH,
    'MultipleWinners'= MultipleWinners,
    'rej_flag_Prev'= rej_flag_Prev ,
    'rej_flag_Curr'= rej_flag_Curr,
    'SelectionLook'= c(),
    'SelectedIndex'= NULL,
    'DropedFlag'=DropedFlag,
    'LastLook'= nLooks,
    'Modify'=F,
    'ModificationLook'=c(),
    'newWeights' =NA,
    'newG' = NA,
    'bdryTab'= NA,
    'HypoMap' = HypoMap,
    'AllocSampleSize'=NA,
    'Stage2AllocSampleSize'=NA,
    'Stage1Obj'=NA,
    'AdaptObj'=NA,
    'allGraphsPrev'=allGraphs,
    'allGraphs'=allGraphs,
    'SubText'=SubText
  )

  look = 1; ContTrial = T
  while(ContTrial) ## Loop for each Interim Analysis
  {
    Prev_mcpObj <- mcpObj #to store inputs before modification

    mcpObj$CurrentLook <- look
    p_raw <- getRawPValues(mcpObj) ## User Input raw p-values

    mcpObj$p_raw <- addNAPvalue(p_raw, GlobalIndexSet)

    if(mcpObj$CurrentLook == 1){

      Stage1Test <- PerformStage1Test(nArms = nArms, nEps = nEps, EpType = EpType, nLooks = nLooks,
                                      nHypothesis = nHypothesis, sigma = sigma,prop.ctr = prop.ctr,
                                      allocRatio = allocRatio,SampleSize = SampleSize,
                                      alpha = alpha, info_frac = info_frac,
                                      typeOfDesign = typeOfDesign,des.type = des.type,
                                      test.type = test.type, Stage1Pvalues = mcpObj$p_raw,
                                      HypoMap = mcpObj$HypoMap, WH = mcpObj$WH)
      cat('Planned Variance Covariance Matrix \n')
      print(Stage1Test$Stage1Obj$Sigma)

      cat("Stage-1 Output Tables \n")
      print(Stage1Test$Stage1Tables)

      mcpObj$rej_flag_Prev <- mcpObj$rej_flag_Curr <- Stage1Test$Stage1Obj$Stage1Analysis$PrimaryHypoTest
      mcpObj$IndexSet <-  paste("H",which(!mcpObj$rej_flag_Curr), sep ='')
      mcpObj$AllocSampleSize <- mcpObj$Stage2AllocSampleSize <- Stage1Test$Stage1Obj$AllocSampleSize
      mcpObj$Stage1Obj <- Stage1Test$Stage1Obj

      if(length(which(mcpObj$rej_flag_Curr==T)) != 0)
      {
        rejected <- paste("H",which(mcpObj$rej_flag_Curr==T), sep ='')
      }else
      {
        rejected <- c()
      }

      mcpObj$WH_Prev <- Stage1Test$Stage1Obj$WH

      if(length(rejected)==0)
      {
        mcpObj$WH <- mcpObj$WH_Prev
      }else
      {
        mcpObj$WH <- mcpObj$WH_Prev[
          apply(mcpObj$WH_Prev[rejected],1, function(x){any(x==1)}) != 1,]
        row.names(mcpObj$WH) <- NULL
      }

      if(plotGraphs) #Plot after Stage-1 analysis
      {
        HypothesisName <- mcpObj$allGraphs$HypothesisName
        HypoIDX <- get_numeric_part(HypothesisName)
        activeStatus <- !unlist(mcpObj$rej_flag_Curr[HypoIDX])
        graphIDX <- which(mcpObj$allGraphs$IntersectIDX == paste(as.integer(activeStatus),collapse = ''))

        if(length(graphIDX)==0)
        {
          nodes <- edges <- NULL
        }else
        {
          nodes <- mcpObj$allGraphs$IntersectionWeights[graphIDX,
                                                        grep('Weight', names(mcpObj$allGraphs$IntersectionWeights))]
          edges <- mcpObj$allGraphs$Edges[[graphIDX]]
        }

        plotGraph(HypothesisName = HypothesisName,
                  w = unlist(nodes),
                  G = edges,
                  activeStatus = activeStatus,
                  Titel = paste('Graph After Stage ', mcpObj$CurrentLook,' analysis'),
                  Text = mcpObj$SubText)
      }

    }else{
      #Stage 2 Analysis
      Stage2Test <- PerformStage2Test(mcpObj=mcpObj, AdaptStage2= AdaptStage2)
      if(AdaptStage2){
        cat('Variance Covariance matrix after adaptation \n')
        print(mcpObj$AdaptObj$Stage2Sigma)
      }
      cat("Stage-2 Output Tables \n")
      print(Stage2Test$Stage2Tables)
      mcpObj$rej_flag_Curr <- Stage2Test$RejStat

      if(plotGraphs) #Plot after Stage-2 analysis
      {
        HypothesisName <- mcpObj$allGraphs$HypothesisName
        HypoIDX <- get_numeric_part(HypothesisName)
        activeStatus <- !unlist(mcpObj$rej_flag_Curr[HypoIDX])
        graphIDX <- which(mcpObj$allGraphs$IntersectIDX == paste(as.integer(activeStatus),collapse = ''))

        if(length(graphIDX)==0)
        {
          nodes <- edges <- NULL
        }else
        {
          nodes <- mcpObj$allGraphs$IntersectionWeights[graphIDX,
                                                        grep('Weight', names(mcpObj$allGraphs$IntersectionWeights))]
          edges <- mcpObj$allGraphs$Edges[[graphIDX]]
        }
        plotGraph(HypothesisName = HypothesisName,
                  w = unlist(nodes),
                  G = edges,
                  activeStatus = activeStatus,
                  Titel = paste('Graph After Stage ', mcpObj$CurrentLook,' analysis'),
                  Text = mcpObj$SubText)
      }
    }

    #Pre-computation for the next look
    # Choice to proceed to next look or start over
    trialTermUserInput <- terminateTrial(mcpObj)

    if(trialTermUserInput == 'y') #proceed to next look
    {
      if(AdaptStage2)
      {
        # Selection for next look
        if(Selection & (look< nLooks))
        {
          mcpObj <- do_Selection(mcpObj)
          if(length(mcpObj$SelectedIndex)!=0) #If the selection set is non empty
          {
            contArms <- getArmsFromHypo(SetH = mcpObj$SelectedIndex, Hypo_map = mcpObj$HypoMap)
            mcpObj$ArmsPresent <- contArms

            if(plotGraphs) #Plot after selection
            {
              HypothesisName <- mcpObj$allGraphs$HypothesisName
              HypoIDX <- get_numeric_part(HypothesisName)
              ActiveIDX <-  get_numeric_part(mcpObj$IndexSet)
              activeStatus <- rep(F,length(HypothesisName))
              activeStatus[which(HypoIDX %in% ActiveIDX)] <- T
              graphIDX <- which(mcpObj$allGraphs$IntersectIDX == paste(as.integer(activeStatus),collapse = ''))

              nodes <- mcpObj$allGraphs$IntersectionWeights[graphIDX,
                                                            grep('Weight', names(mcpObj$allGraphs$IntersectionWeights))]
              edges <- mcpObj$allGraphs$Edges[[graphIDX]]
              plotGraph(HypothesisName = HypothesisName,
                        w = unlist(nodes),
                        G = edges,
                        activeStatus = activeStatus,
                        Titel = paste('Graph After Selection'),
                        Text = mcpObj$SubText)

            }
          }else #If the selection set is empty
          {
            mcpObj$ContTrial = F
          }
        }

        # Stage-2 Samples
        if(ModifySamples)
        {
          Stage2ModSS <- do_ModifyStage2Sample(allocRatio=mcpObj$allocRatio,
                                               ArmsPresent = mcpObj$ArmsPresent,
                                               AllocSampleSize = mcpObj$AllocSampleSize)

          mcpObj$Stage2AllocSampleSize <- Stage2ModSS$newAllocSampleSize
          mcpObj$Stage2allocRatio <- Stage2ModSS$newallocRatio
        }

        #Modify the weights and testing strategy
        if(UpdateStrategy & (look< nLooks) & (length(mcpObj$IndexSet)>1))
        {
          mcpObj <- do_modifyStrategy(mcpObj,showExistingStrategy = F)
        }

        #Modify the Stage2 boundaries
        AdaptResults <- adaptBdryCER(mcpObj)
        mcpObj$AdaptObj <- AdaptResults
      }

      look <- look + 1

    }else if(trialTermUserInput == 'n') #terminate the trial
    {
      break

    }else if(trialTermUserInput == 's') #Start over from the last look inputs
    {
      mcpObj <- Prev_mcpObj
    }
    #------------------End of the while loop------------------------#
  }
}
