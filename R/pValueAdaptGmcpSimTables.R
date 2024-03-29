
# The file contains supporting functions(Detailed Output Tables) for Adaptive GMCP Simulation Function AdaptGmcp_Simulation/adaptGMCP_SIM(.) ----
## Author: Ajoy.M

#-------------- -
#Compute Overall Power Table from all the simulation
#-------------- -
SimPowers <- function(nSimulation,PowerTab)
{
  values <- as.numeric(apply(PowerTab[,-1], 2, function(x){sum(x)/nSimulation}))
  PowerTable <- data.frame('Overall_Powers'=c('Global Power','Conjunctive Power','Disjunctive Power','FWER'),
                           'Values'= values)
  PowerTable
}

#------------- -
#Count contribution to different powers from each simulations
#------------- -
CountPower <- function(simID, SummaryStatFile,TrueNull)
{
  rejMat <- subset(SummaryStatFile,SimID==simID)
  rejMat <- rejMat[,grep('RejStatus', names(rejMat))]
  rej.final <- apply(rejMat, 2, function(col) any(col, na.rm = T))
  data.frame(
    'simID' = simID,
    'nG' = as.integer(any(rej.final)),
    'nC' = as.integer(
      ifelse(length(rej.final[!TrueNull])==0,0,all(rej.final[!TrueNull]))
      ),
    'nD' = as.integer(
      ifelse(length(rej.final[!TrueNull])==0,0,any(rej.final[!TrueNull]))
      ),
    'nF' = as.integer(
      ifelse(length(rej.final[TrueNull])==0,0,any(rej.final[TrueNull]))
      )
  )
}

#------------ -
#Identify True Null based of the response generation inputs
#------------ -
checkTrueNull <- function(gmcpSimObj,index_map)
{
  pairsHyp <- strsplit(index_map$Names,split = '-')

  delta <- unlist(lapply(pairsHyp, function(x){
    abs(gmcpSimObj$Arms.Mean[x[1]] - gmcpSimObj$Arms.Mean[x[2]])
  }))
  index_map$TrueNull <- (delta<1E-6)
  index_map

}

#------------ -
#Identify True Null based of the response generation inputs
#------------ -
checkTrueNull2 <- function(HypoMap,Arms.Mean)
{
  delta <- unlist(lapply(1:nrow(HypoMap), function(i){
    abs(Arms.Mean[[HypoMap$Groups[i]]][HypoMap$Treatment[i]] -
          Arms.Mean[[HypoMap$Groups[i]]][HypoMap$Control[i]])
  }))
  TrueNull <- (delta<1E-6)
  TrueNull
}

#------------ -
#Identify True Null based of the response generation inputs
#------------ -
checkTrueNull3 <- function(HypoMap,Arms.Mean,Arms.Prop)
{
  delta <- unlist(lapply(1:nrow(HypoMap), function(i){
    if(HypoMap$EpType[i] == 'Continuous'){

      abs(Arms.Mean[[HypoMap$Groups[i]]][HypoMap$Treatment[i]] -
            Arms.Mean[[HypoMap$Groups[i]]][HypoMap$Control[i]])

    }else if(HypoMap$EpType[i] == 'Binary'){

      abs(Arms.Prop[[HypoMap$Groups[i]]][HypoMap$Treatment[i]] -
            Arms.Prop[[HypoMap$Groups[i]]][HypoMap$Control[i]])

    }

  }))
  TrueNull <- (delta<1E-6)
  TrueNull
}

#------------- -
#Count contribution to different powers from each simulations
#------------- -
CountEfficacy <- function(simID, SummaryStatFile)
{
  rejMat <- subset(SummaryStatFile,SimID==simID)
  rejMat <- rejMat[,grep('RejStatus', names(rejMat))]

  rej.final <- apply(rejMat, 2, function(col) any(col, na.rm = T))
  n <- length(rej.final)
  if(n==0) rej.final <- rep(F, n)
  interHypo <- genCombs(n)

  idx <- lapply(1:nrow(interHypo), function(x)paste(interHypo[x,], collapse =''))
  rej_idx <- paste(as.integer(rej.final), collapse ='')

  col_name <- lapply(1:nrow(interHypo),function(x){
    paste(paste('H',which(interHypo[x,] == 1),sep = ''), collapse = ',')
  })

  eff <- rep(0, length(col_name))
  eff[which(idx == rej_idx)] <- 1

  eff_count <- data.frame(matrix(c(simID,eff),nrow = 1))
  colnames(eff_count) <- c('simID', col_name)
  eff_count
}

genCombs <- function(n) {
  combs <- rep(0, n)
  combs <- data.frame(do.call(rbind,lapply(0:n, function(i)
                     t(apply(combn(1:n,i), 2, function(k) {combs[k]=1;combs})))))
  combs <- combs[-1,]
  rownames(combs) <-  NULL
  return(combs)
}


