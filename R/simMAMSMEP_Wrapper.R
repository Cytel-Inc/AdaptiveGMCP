#' Wrapper function that takes all inputs through a single R dataframe
#' @param InputDF R Dataframe: This is the csv/excel input data in the R dataframe format
#' @export
simMAMSMEP_Wrapper <- function(InputDF){
  # Update the dataframe column names in the following mapping in case
  # the names in the input csv/excel changes
  library(dplyr)
  library(tidyr)

  lOut <- list()
  for (nModelNum in 1:nrow(InputDF)){

    out <- tryCatch({run1TestCase(InputDF = InputDF[nModelNum, ])},
                    error = function(err){
                      paste0("Model ", nModelNum, " execution failed.")
                    })
    errorLog <- c(!grepl('Invalid', out[[1]]), !is.character(out), !is.null(out))
    if(!any(errorLog == F)){
      # extract the power table from the output
      dfOverall_Powers_long <- out$Overall_Powers_df
      # convert the power table to wide format and add serial number column
      dfOverall_Powers_wide <- dfOverall_Powers_long %>%
        pivot_wider(names_from = Overall_Powers, values_from = Values) %>%
        mutate(Sno = nModelNum) %>%
        relocate(Sno, .before = everything())
      #Model ID and Seed Number to reproduce
      mInfo <- data.frame('ModelID'=InputDF[nModelNum,'ModelID'], 'Seed'=out$Seed)

      #Output Table
      OutTab <- cbind(mInfo, data.frame(dfOverall_Powers_wide[,-1]))
      # add each iterations power table to a list
      lOut[[nModelNum]] <- OutTab

      passedTxt <- paste0("Model ", nModelNum, " execution completed successfully.")
      cat('\n',passedTxt,'\n')
      print(paste0("Power table for ", nModelNum, ":"))
      print(OutTab)


    }else if(grepl('Invalid', out[[1]])){
      failTxt <- paste0("Model ", nModelNum, " execution failed.")
      cat('\n',failTxt,'\n')
      cat('\n Details \n')
      print(unlist(out))
    }else{
      print(out)
    }
  }
  # rbind power tables for each iteration to produce a single table
  dfOut <- do.call(rbind, lOut)
  return(dfOut)
}


run1TestCase <- function(InputDF){
  # mapping to link simMAMSMEP function arguments with csv columns
  Method <- InputDF$Method
  SampleSize <- InputDF$SampleSize
  alpha <- InputDF$alpha
  TestStat <- InputDF$TestStat
  FWERControl <- InputDF$FWERControl
  nArms <- InputDF$nArms
  nEps <- InputDF$nEps
  Arms.Mean <- eval(parse(text = InputDF$Arms.Mean))
  Arms.std.dev <- eval(parse(text = InputDF$Arms.std.dev))
  Arms.alloc.ratio <- eval(parse(text = InputDF$Arms.alloc.ratio))
  EP.Corr <- eval(parse(text = InputDF$EP.Corr))
  WI <- eval(parse(text = InputDF$WI))
  G <- eval(parse(text = InputDF$G))
  test.type <- InputDF$test.type
  info_frac <- eval(parse(text = InputDF$info_frac))
  typeOfDesign <- InputDF$typeOfDesign
  MultipleWinners <- InputDF$MultipleWinners
  Selection <- InputDF$Selection
  SelectionLook <- InputDF$SelectionLook
  SelectEndPoint <- InputDF$SelectEndPoint
  SelectionScale <- InputDF$SelectionScale
  SelectionCriterion <- InputDF$SelectionCriterion
  SelectionParmeter <- InputDF$SelectionParmeter
  KeepAssosiatedEps <- InputDF$KeepAssosiatedEps
  ImplicitSSR <- InputDF$ImplicitSSR
  nSimulation <- InputDF$nSimulation
  Seed <- InputDF$Seed
  SummaryStat <- InputDF$SummaryStat
  plotGraphs <- InputDF$plotGraphs
  Parallel <- InputDF$Parallel
  # put the following code in try catch so the loop continues even if one iteration fails
  out <- simMAMSMEP(alpha = alpha, SampleSize = SampleSize,nArms = nArms,nEps = nEps, TestStat=TestStat, FWERControl = FWERControl,

                    Arms.Mean = Arms.Mean, Arms.std.dev = Arms.std.dev, Arms.alloc.ratio = Arms.alloc.ratio,

                    EP.Corr = EP.Corr,WI = WI, G = G, test.type = test.type,info_frac = info_frac,

                    typeOfDesign = typeOfDesign, MultipleWinners = MultipleWinners,

                    Selection = Selection,SelectionLook = SelectionLook,SelectEndPoint = SelectEndPoint,SelectionScale = SelectionScale,

                    SelectionCriterion = SelectionCriterion, SelectionParmeter = SelectionParmeter, KeepAssosiatedEps = KeepAssosiatedEps,

                    ImplicitSSR = ImplicitSSR, nSimulation = nSimulation, Seed = Seed, SummaryStat = SummaryStat,

                    Method = Method, plotGraphs = plotGraphs, Parallel = Parallel)
  out
}


#' Create table and plots of given format
#' @param dfOut : output from simMAMSMEP_Wrapper
#' @export
genPowerTablePlots <- function(PowerType, dfOut, TableTemDF){
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  library(tidyr)

  p <- sapply(TableTemDF$ModelID, function(mID){
    p <- dfOut[dfOut$ModelID == mID, PowerType]
    if(length(p)==0){
      p <- NA
    }
    p
  })
  tab1 <- TableTemDF[,- which(names(TableTemDF)=='ModelID')]
  tab1[PowerType] <- p
  scenarios <- unique(tab1$Level1)
  plots <- list()

  # for(s in 1:length(scenarios)){
  #   tab2 <- tab1[tab1$Level1 == scenarios[s], ]
  #   tab2[PowerType][is.na(tab2[PowerType])] = 0
  #   plots[[s]] <- ggplot(tab2, aes(x = Level2)) +
  #     geom_bar(aes(y = Global.Power, fill = Method), stat = "identity", position = "dodge") +
  #     ggtitle(scenarios[s])+
  #     theme(plot.title = element_text(hjust = 0.5),
  #           axis.title.x = element_blank())
  #
  # }
  # plt <- grid.arrange(grobs=plots,ncol=2)
  tab3 <- reshape(data = tab1, idvar = c('Level1','Level2'),
          v.names = PowerType,timevar = 'Method',direction = 'wide')
  tab3$Difference <- tab3[,4]- tab3[,3]
  colnames(tab3) <- c('Scenario',
                      'Treatment Selection Rule',
                      'Stagewise MAMS',
                      'cumulative MAMS',
                      'Difference')
  SelectionRuleOrder <- c("Conservative", "Normal", "Aggressive", "Ultra")
  tab3$`Treatment Selection Rule` <- factor(tab3$`Treatment Selection Rule`, levels = SelectionRuleOrder)
  tabLong <- tab3 %>%
    select(!Difference) %>%
    pivot_longer(cols = !c('Scenario', 'Treatment Selection Rule'),
                 names_to = 'MAMS',
                 values_to = 'value',
                 values_drop_na = TRUE)

 list('TableWide' = tab3, 'TableLong' = tabLong)
}
