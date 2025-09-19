# Batch_gMCP_FS_Sim_Bin.R
# File for executing a batch of gMCP simulation scenarios

library(tidyverse)
library(gMCP)
# source("R/gMCP_calcPower_wrapper.R", keep.source = TRUE)

dfInput <- read_csv("internalData/BatchInp_gMCP_FS_Sim_Bin.csv")
sOutFilePrefix <- "Output"
sOutPath <- "internalData/"

nModelsToRun <- dfInput$ModelID

# TRIAL RUN - START >>>>>>>>>>>>
# To do a trial run, uncomment this block so that the tests are run with a
# small number of simulations rather than the number specified in the input
# file.
dfInput$nSimulation <- 10 # 10000
# TRIAL RUN - OVER >>>>>>>>>>>>

tStartTime <- Sys.time()

dfOutput <- gMCP_CalcPowerWrapper(InputDF = dfInput %>%
                                     filter(ModelID %in% nModelsToRun))

# dfOutput <- lOut$gMCPOut
# gMCPInp <- lOut$gMCPInp

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

