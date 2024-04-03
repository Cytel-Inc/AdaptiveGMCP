
#Load the package
library(AdaptGMCP)

# Input Csv containing all scenarios
InputDF <- read.csv('internalData/TestCases/InputScenariosBinaryNormal.csv',
                    fileEncoding = "UTF-8-BOM")

#Set Parallel as TRUE(always set as TRUE)
InputDF$Parallel <- TRUE

#Specify the modelIDs to run
#If all models need to be executed then set it as  modelsToRun <- InputDF$ModelID
modelsToRun <- InputDF$ModelID

start.time <- Sys.time()
outDF <- simMAMSMEP_Wrapper(InputDF = InputDF %>% filter(ModelID %in% modelsToRun))
elapsed.time <- Sys.time()-start.time

#Print Output
outDF

# #Save Output
# timeNow <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
# write.csv(outDF, paste0("Outputs/outDF_",timeNow,".csv"), row.names = F)
#
# #Execution details
# SysInfo <- Sys.info()
# cat("Execution performed in ", SysInfo['nodename'],'\n',
#     "Execution time ", elapsed.time)
#
