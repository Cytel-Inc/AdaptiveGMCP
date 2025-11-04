# File: Batch_2Lk_AdaptGMCP_Sim_Bin.R
# This file contains code for executing a batch of simulations for 2 look
# GMCP designs with binary endpoints.

library(tidyverse)

# We will use the function simMAMSMEP_Wrapper() for executing a batch.

# dfInput <- read_csv("internalData/InputScenarios_2ep5arm.csv")
dfInput <- read_csv("internalData/CER_Inp_1ep5arms.csv")
sOutFilePrefix <- "AGMCP_2Lk_Out"
sOutPath <- "internalData/"

nModelsToRun <- c(1, 2, 3) # dfInput$ModelID

# TRIAL RUN - START >>>>>>>>>>>>
# To do a trial run, uncomment this block so that the tests are run with a
# small number of simulations rather than the number specified in the input
# file.
dfInput$nSimulation <- 10 #
dfInput$Parallel <- FALSE
# dfInput$test.type <- "Parametric"
# TRIAL RUN - OVER >>>>>>>>>>>>

tStartTime <- Sys.time()

dfOutput <- simMAMSMEP_Wrapper(InputDF = dfInput %>%
                                 filter(ModelID %in% nModelsToRun))

tElapTime <- Sys.time() - tStartTime

dfOutput

#Save Output
sTimeNow <- format(Sys.time(), "%d%h%y-%H_%M")
sOutPath1 <- paste0(sOutPath, sOutFilePrefix, "_", sTimeNow, ".csv")
# sOutPath1 <- paste0(sOutPath, "Output_FS_GMCP_Sim_Bin", ".csv")
write.csv(dfOutput, sOutPath1, row.names = F)

#Execution details
SysInfo <- Sys.info()
cat("Execution performed on ", SysInfo['nodename'],"\n",
    "Execution time", tElapTime)
