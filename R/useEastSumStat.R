# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

#To use East Summary Statistics File(only applicable for one look designs)
useEastSumStat <- function(SimID, EastSumStat){
  #The following is specific to East two sample mixed type designs(may not work for any other designs)
  EastSumStatSim <- EastSumStat[EastSumStat$Sim_Index==SimID ,]
  LookID <- 1 #always 1 for FSD

  Delta <- EastSumStatSim[,grep('Delta',names(EastSumStatSim))]
  names_Delta <- paste('Delta', 1:length(Delta), sep = '')

  StdError <- EastSumStatSim[,grep('StdError',names(EastSumStatSim))]
  names_StdError <- paste('StdError', 1:length(StdError), sep = '')

  TestStat <- EastSumStatSim[,grep('Test_Stat',names(EastSumStatSim))]
  names_TestStat <- paste('TestStat', 1:length(TestStat), sep = '')

  RawPvalues <- EastSumStatSim[,grep('pval_1_Endpoint',names(EastSumStatSim))]
  names_RawPvalues <- paste('RawPvalues', 1:length(RawPvalues), sep = '')

  SumStat <- data.frame(matrix(c(SimID,LookID,
                                 Delta,StdError,
                                 TestStat,RawPvalues),nrow = 1))
  colnames(SumStat) <- c("SimID","LookID",
                         names_Delta,names_StdError,
                         names_TestStat,names_RawPvalues)
  SumStat
}



