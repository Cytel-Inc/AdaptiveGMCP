
#Uncomment the following code to run simulation;
#short cut: 1) select all the following lines 2) ctrl+shift+c]


#Load the package
library(AdaptGMCP)

###########-----------Inputs Need to be Modified for each run------#########
# Input Csv containing all scenarios(Change the path before reading the data)
InputDF <- read.csv('internalData/TestCases/InputScenariosBinaryNormal.csv',
                    fileEncoding = "UTF-8-BOM")

# Output csv file name(without the .csv extension)
# Format <OutFileName>_<execution date and time>.csv
OutFileName <- "Output"

#Output file directory
OutPath <- "C:\\Users\\neeraj\\OneDrive - Cytel\\Documents\\CIC\\Code\\OutPuts_From_Wrapper_Fun\\"

##################################################################################
###########-----------Inputs Need to be Modified for customized run------#########

# Set Parallel as TRUE(always set as TRUE)
InputDF$Parallel <- T

# Specify the modelIDs to run
# If all models need to be executed then set it as  modelsToRun <- InputDF$ModelID
modelsToRun <- InputDF$ModelID

##################################################################################
##########-----------Execution, No Modification is required-----------############

#Start Time
start.time <- Sys.time()

#Run Simulation Batches

outDF <- simMAMSMEP_Wrapper(InputDF = InputDF %>% filter(ModelID %in% modelsToRun))

#Elapsed Time
elapsed.time <- Sys.time()-start.time

#Print Output
outDF

#Save Output
timeNow <- format(Sys.time(), "%d%h%y-%H_%M")
OutPath1 <- paste0(OutPath,OutFileName,"_",timeNow,".csv")
write.csv(outDF,OutPath1)

#Execution details
SysInfo <- Sys.info()
cat("Execution performed in ", SysInfo['nodename'],'\n',
    "Execution time ", elapsed.time)
###################################################################################

