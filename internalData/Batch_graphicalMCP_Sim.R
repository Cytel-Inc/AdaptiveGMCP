# File: Batch_graphicalMCP_Sim.R

source("internalData/graphicalMCP_pow_calc_wrapper.R")

# Read graph input data
dfInput <- read_csv("internalData/BatchInput_AGMCP_Sim_Bin_3.csv")
sOutFilePrefix <- "grMCP_Out"
sOutPath <- "internalData/"

# Filtering out Sidak test as graphicalMCP does not appear to support it
dfInput <- dfInput %>% filter(test.type != "Sidak")

nModelsToRun <- dfInput$ModelID # c(1, 2, 3) # c(110, 112, 113, 115) #

# TRIAL RUN - START >>>>>>>>>>>>
# To do a trial run, uncomment this block so that the tests are run with a
# small number of simulations rather than the number specified in the input
# file.
dfInput$nSimulation <- 100 # 10000 # 10
dfInput$Parallel <- FALSE
# TRIAL RUN - OVER >>>>>>>>>>>>


tStartTime <- Sys.time()

dfOutput <- graphicalMCP_Wrapper(InputDF = dfInput %>%
                                 filter(ModelID %in% nModelsToRun))

tElapTime <- Sys.time() - tStartTime

dfOutput

#Save Output
sTimeNow <- format(Sys.time(), "%d%h%y-%H_%M")
sOutPath1 <- paste0(sOutPath, sOutFilePrefix, "_", sTimeNow, ".csv")
write.csv(dfOutput, sOutPath1, row.names = F)

#Execution details
SysInfo <- Sys.info()
cat("Execution performed on ", SysInfo['nodename'],"\n",
    "Execution time", tElapTime)
