# File: Batch_2Lk_AdaptGMCP_Sim_Bin.R
# This file contains code for executing a batch of simulations for 2 look
# GMCP designs with binary endpoints.

# library(tidyverse)

nSim <- 10 # 1000 # 10000 # 50000
nSim2 <- 5 # 100 # 1000 # 1000

out <- simMAMSMEP(
  Method = "CER", alpha = 0.025, SampleSize = 162, nArms = 3, nEps = 1,
  lEpType=list('EP1' = 'Binary'), TestStatBin = "UnPooled", FWERControl = "CombinationTest",
  Arms.Prop = list('EP1' = c(0.1, 0.1, 0.1)),
  Arms.alloc.ratio = c(1,1, 1), EP.Corr = matrix(1), WI = c(rep(1/2,2)),
  G = rbind(H1=c(0,1), H2=c(1,0)), test.type = "Partly-Parametric",
  info_frac = c(0.75,1), typeOfDesign = "asOF",
  MultipleWinners = T, Selection = F, SelectionLook = NA, SelectEndPoint = NA,
  SelectionScale = NA, SelectionCriterion = NA, SelectionParmeter = NA,
  KeepAssosiatedEps = NA, ImplicitSSR = "None",
  nSimulation = nSim, nSimulation_Stage2 = nSim2, Seed = 1234, SummaryStat = T,
  plotGraphs = F, Parallel = T # F
)

print(out)


### HELP EXAMPLE ##############################################
    # Method <- 'CER'
    # SampleSize <- 400
    # alpha <- 0.025
    # nArms <- 4
    # nEps  <- 2
    # # nArms <- 3
    # # nEps  <- 3
    # # EpType <- list("EP1" = "Continuous",
    # #                "EP2" = "Continuous")
    # # EpType <- list("EP1" = "Binary",
    # #                "EP2" = "Binary")
    # EpType <- list("EP1" = "Continuous",
    #                "EP2" = "Binary")
    # # EpType <- list("EP1" = "Binary",
    # #                "EP2" = "Continuous",
    # #                "EP3" = "Binary")
    # TestStatCon <- "t-equal"
    # TestStatBin <- "Pooled" # "UnPooled"
    # FWERControl <- "CombinationTest" # None"
    # # Arms.Mean <- list('EP1' = c(0, 0.4, 0.4, 0.4),
    # #                   'EP2' = c(0, 0.4, 0.4, 0.4))
    # # Arms.Mean <- list('EP1' = NA,
    # #                   'EP2' = c(0, 0.4, 0.4),
    # #                   'EP3' = NA)
    # # Arms.Mean <- list('EP1' = NA, 'EP2' = NA)
    # Arms.Mean <- list('EP1' = c(0, 0.4, 0.4, 0.4), 'EP2' = NA)
    # # Arms.std.dev <- list('EP1' = c(1.1, 1.2, 1.3, 1.4),
    # #                      'EP2' = c(1.1, 1.2, 1.3, 1.4))
    # # Arms.std.dev <- list('EP1' = NA,
    # #                      'EP2' = c(1.1, 1.2, 1.3),
    # #                      'EP3' = NA)
    # # Arms.std.dev <- list('EP1' = NA, 'EP2' = NA)
    # Arms.std.dev <- list('EP1' = c(1.1, 1.2, 1.3, 1.4), 'EP2' = NA)
    #
    # # Arms.Prop <- list(EP1 = NA, EP2 = NA)
    # # Arms.Prop <- list(EP1 = c(0.2, 0.35, 0.45), EP2 = NA, EP3 = c(0.1, 0.1, 0.1))
    # # Arms.Prop <- list(EP1 = c(0.2, 0.35, 0.45, 0.2), EP3 = c(0.1, 0.1, 0.1, 0.1))
    # Arms.Prop <- list(EP1 = NA, EP2 = c(0.2, 0.35, 0.45, 0.2))
    #
    # CommonStdDev <- FALSE
    # Arms.alloc.ratio <- c(1, 1, 1, 1)
    # # Arms.alloc.ratio <- c(1, 1, 1)
    # EP.Corr <- matrix(c(1, 0.5,
    #                     0.5, 1),
    #                   nrow = nEps)
    # # EP.Corr <- matrix(c(1, 0.5, 0.5,
    # #                     0.5, 1, 0.5,
    # #                     0.5, 0.5, 1),
    # #                   nrow = nEps)
    # WI <-  c(1/3, 1/3, 1/3, 0, 0, 0)
    # G <- matrix(c(0,0,0,     1,0,0,
    #               0,0,0,     0,1,0,
    #               0,0,0,     0,0,1,
    #               0,1/2,1/2, 0,0,0,
    #               1/2,0,1/2, 0,0,0,
    #               1/2,1/2,0, 0,0,0),
    #             nrow = nEps*(nArms-1), byrow = TRUE)
    # test.type <- 'Partly-Parametric'
    # info_frac <-  c(1/2,1)
    # typeOfDesign <- "asOF"
    # MultipleWinners <- TRUE
    # Selection <- FALSE # TRUE
    # SelectionLook <- 1
    # SelectEndPoint <- 1
    # SelectionScale <- 'teststat'
    # SelectionCriterion <- 'threshold'
    # SelectionParmeter <- 0.6745
    # KeepAssosiatedEps <- TRUE
    # ImplicitSSR <- 'Selection'
    # nSimulation <- 10
    # nSimulation_Stage2 <- 10
    # Seed <- 100
    # SummaryStat <- FALSE
    # plotGraphs <- FALSE
    # Parallel <- FALSE
    # out <- simMAMSMEP(
    #   alpha = alpha, SampleSize = SampleSize, nArms = nArms, nEps = nEps,lEpType=EpType,
    #   TestStatCon = TestStatCon, TestStatBin = TestStatBin, FWERControl = FWERControl,
    #   Arms.Mean = Arms.Mean, Arms.std.dev = Arms.std.dev, CommonStdDev = CommonStdDev,
    #   Arms.Prop = Arms.Prop, Arms.alloc.ratio = Arms.alloc.ratio,
    #   EP.Corr = EP.Corr, WI = WI, G = G, test.type = test.type, info_frac = info_frac,
    #   typeOfDesign = typeOfDesign, MultipleWinners = MultipleWinners,
    #   Selection = Selection, SelectionLook = SelectionLook, SelectEndPoint = SelectEndPoint, SelectionScale = SelectionScale,
    #   SelectionCriterion = SelectionCriterion, SelectionParmeter = SelectionParmeter, KeepAssosiatedEps = KeepAssosiatedEps,
    #   ImplicitSSR = ImplicitSSR, nSimulation = nSimulation, nSimulation_Stage2 = nSimulation_Stage2, Seed = Seed, SummaryStat = SummaryStat,
    #   Method = Method, plotGraphs = plotGraphs, Parallel = Parallel
    # )
    #
    # outPower <- out$Overall_Powers
    # outPower
# HELP EXAMPLE OVER
###############################################################

### 1 CONTINUOUS EP TEST ##########################
## WITH TREATMENT SELECTION
# out <- simMAMSMEP(
#   Method = "CER", alpha = 0.025, SampleSize = 500, nArms = 5, nEps = 1,
#   lEpType=list('EP1' = 'Continuous'), TestStatCon = "t-equal", FWERControl = "None",
#   Arms.Mean = list('EP1' = c(0, 0, 0, 0, 0.25)),
#   Arms.std.dev = list('EP1' = c(1, 1, 1, 1, 1)), CommonStdDev = F,
#   Arms.alloc.ratio = c(1,0.5,0.5,0.5,0.5), EP.Corr = matrix(1), WI = c(rep(1/4,4)),
#   G = rbind(H1=c(0,1/3,1/3,1/3), H2=c(1/3,0,1/3,1/3), H3=c(1/3,1/3,0,1/3), H4=c(1/3,1/3,1/3,0)),
#   test.type = "Parametric", info_frac = c(1/2,1), typeOfDesign = "asOF",
#   MultipleWinners = T, Selection = T, SelectionLook = 1, SelectEndPoint = 1,
#   SelectionScale = "pvalue", SelectionCriterion = "threshold",
#   SelectionParmeter = 0.75, KeepAssosiatedEps = T, ImplicitSSR = "Selection",
#   nSimulation = 3, nSimulation_Stage2 = 100, Seed = 1234, SummaryStat = T,
#   plotGraphs = F, Parallel = F
# )

## WITHOUT TREATMENT SELECTION
# out <- simMAMSMEP(
#   Method = "CER", alpha = 0.025, SampleSize = 500, nArms = 5, nEps = 1,
#   lEpType=list('EP1' = 'Continuous'), TestStatCon = "t-equal", FWERControl = "None",
#   Arms.Mean = list('EP1' = c(0, 0, 0, 0, 0.25)),
#   Arms.std.dev = list('EP1' = c(1, 1, 1, 1, 1)), CommonStdDev = F,
#   Arms.alloc.ratio = c(1,0.5,0.5,0.5,0.5), EP.Corr = matrix(1), WI = c(rep(1/4,4)),
#   G = rbind(H1=c(0,1/3,1/3,1/3), H2=c(1/3,0,1/3,1/3), H3=c(1/3,1/3,0,1/3), H4=c(1/3,1/3,1/3,0)),
#   test.type = "Parametric", info_frac = c(1/2,1), typeOfDesign = "asOF",
#   MultipleWinners = T, Selection = F, SelectionLook = NA, SelectEndPoint = NA,
#   SelectionScale = NA, SelectionCriterion = NA, SelectionParmeter = NA,
#   KeepAssosiatedEps = NA, ImplicitSSR = "None",
#   nSimulation = 3, nSimulation_Stage2 = 100, Seed = 1234, SummaryStat = T,
#   plotGraphs = F, Parallel = F
# )
#
# out
#
# ###################################################
# # We will use the function simMAMSMEP_Wrapper() for executing a batch.
#
# # dfInput <- read_csv("internalData/CER_Inp_1ep5arms - Continuous.csv")
# # dfInput <- read_csv("internalData/InputScenarios_2ep5arm.csv")
# # dfInput <- read_csv("internalData/CER_Inp_1ep5arms.csv")
# dfInput <- read_csv("internalData/Inp_CER_Bin_1ep3arms.csv")
# sOutFilePrefix <- "Out_CER_Bin_1ep3arms"
# sOutPath <- "internalData/"
#
# nModelsToRun <- dfInput$ModelID # 36 # 2 # 1 #
#
# # TRIAL RUN - START >>>>>>>>>>>>
# # To do a trial run, uncomment this block so that the tests are run with a
# # small number of simulations rather than the number specified in the input
# # file.
# # dfInput$nSimulation <- 25 #
# # dfInput$nSimulation_Stage2 <- 1000
# # dfInput$Parallel <- FALSE
# # dfInput$test.type <- "Parametric"
# df <- dfInput %>% filter(ModelID %in% nModelsToRun)
# print(paste0("Method = ", df$Method[1])) # CER / CombPValue
# print(paste0("nSimulation = ", df$nSimulation[1]))
# print(paste0("nSimulation_Stage2 = ", df$nSimulation_Stage2[1]))
# print(paste0("Parallel = ", df$Parallel[1]))
# print(paste0("lEpType = ", df$lEpType[1]))
# print(paste0("test.type = ", df$test.type[1])) # Non-Parametric / Bonf
# print(paste0("Seed = ", df$Seed[1]))
# # TRIAL RUN - OVER >>>>>>>>>>>>
#
# tStartTime <- Sys.time()
#
# dfOutput <- simMAMSMEP_Wrapper(InputDF = dfInput %>%
#                                  filter(ModelID %in% nModelsToRun))
#
# tElapTime <- Sys.time() - tStartTime
#
# dfOutput
#
# #Save Output
# sTimeNow <- format(Sys.time(), "%d%h%y-%H_%M")
# sOutPath1 <- paste0(sOutPath, sOutFilePrefix, "_", sTimeNow, ".csv")
# # sOutPath1 <- paste0(sOutPath, "Output_FS_GMCP_Sim_Bin", ".csv")
# write.csv(dfOutput, sOutPath1, row.names = F)
#
# #Execution details
# SysInfo <- Sys.info()
# cat("Execution performed on ", SysInfo['nodename'],"\n",
#     "Execution time", tElapTime)
