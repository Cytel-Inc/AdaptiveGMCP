

#Function to compute the Stage-wise planned boundaries for all intersection
#' @export
getPlanNonParmBdry <- function(nLooks, sig_level, info_frac, typeOfDesign = 'asOF')
{
  if(sig_level != 0)
  {
    des <- rpact::getDesignGroupSequential(kMax = nLooks, alpha = sig_level,
                                           informationRates =info_frac,
                                           typeOfDesign = typeOfDesign)
    if(nLooks == 2)
    {
      return(list('Stage1Bdry' = des$stageLevels[1], 'Stage2Bdry' = des$stageLevels[2]))
    }else if(nLooks == 1)
    {
      return(list('Stage1Bdry' = des$stageLevels[1], 'Stage2Bdry' = 0))

    }else
    {
      print('Error in getPlanNonParmBdry: Boundary Computation is not available for Stages > 2')
    }
  }else
  {
    return(list('Stage1Bdry' = 0, 'Stage2Bdry' = 0))
  }
}

#Component-wise PCER
getPCER <- function(a2, p1, ss1,ss2)
{
  #Make Seed restricted to local##
  old <- .Random.seed
  on.exit( { .Random.seed <<- old } )
  set.seed(200295)
  #################################
  r <- sqrt(ss1/ss2)
  1-pnorm(qnorm(1-a2),mean = r*qnorm(1-p1), sd = sqrt(1-r^2))
}

#Function to compute the stage-2 boundary crossing probability
exitProbStage2Nparam <- function(aj2, aj1,  ss1, ss2) #e: sig level, a1: stage-1 boundary
{
  #Set local seed##
  set.seed(200295)
  #################################
  upper <- c(qnorm(1-aj1), Inf)
  lower <- c(-Inf, qnorm(1-aj2))
  r <- sqrt(ss1/ss2)
  sigma <- matrix(c(1,r, r,1), nrow = 2)
  prob <- mvtnorm::pmvnorm(lower = lower, upper = upper, sigma = sigma)[1]
  return(prob+aj1)
}

#Optimization function to get stage-2 boundary
getBdryStage2Nparam <- function(ej, aj1,  ss1, ss2)
{
  if(ss1 == 0 || ss2 == 0) stop('Error: ss1 == 0 || ss2 == 0 is not true | function: getBdryStage2Nparam')
  minbdry <- 0; maxbdry <- 1 #Set interval
  bdry2NP <-function(x)
  {
    if(ej==1){ #when the the threshold is 1 the exit prob can be assumed to be 1
      extProb <- 1
    }else
    {
      extProb <- exitProbStage2Nparam(x, aj1,  ss1, ss2)
    }
    #cat('aj2 :',x,'extProb:',extProb,'\n')
    extProb - ej
  }

  if((ej-aj1)>0 || ej != 0) #compute boundary when the exit prob >0
  {
    uniroot(f = bdry2NP, interval = c(minbdry, maxbdry), tol = 1E-16)$root

  }else
  {
    0
  }

}

#Function to compute the Adj Boundary based on modified weights
getStage2CondNParamBdry <- function(a1,p1,v,BJ,SS1, SS2)
{
  #Function to compute adjusted PCER for arbtrary g
  getAdjPCER <- function(g,a1,p1, v,SS1, SS2)
  {
    e <- v*g
    a2adj <- unlist(lapply(1:length(e),
                           function(x){getBdryStage2Nparam(ej = e[x], aj1 = a1[x],
                                                     ss1 = SS1[x], ss2 = SS2[x])}
    ))
    A_adj <- unlist(lapply(1:length(a2adj), function(x){
      getPCER(a2 = a2adj[x], p1 = p1[x], ss1 = SS1[x], ss2 = SS2[x])
    }))
    A_adj
  }

  #Function to search for optimum g
  OptimGamma <- function(x)
  {
    if(x==0)
    {
      modBJ <- 0
    }else if( x == 1)
    {
      modBJ <- 1
    }else
    {
      modBJ <- sum(getAdjPCER(g = x,a1 = a1,p1 = p1, v = v, SS1=SS1, SS2=SS2))
    }
    #cat('gammaJ : ' , x, 'BJ : ',modBJ,'\n')
    modBJ-BJ
  }

  #Optimization
  gOpt <- uniroot(OptimGamma,interval = c(0,1),tol = 1E-16)$root

  PCER_adj <- getAdjPCER(g = gOpt,a1 = a1, p1 = p1, v=v, SS1=SS1, SS2=SS2)

  Stage2AdjBdry <- unlist(lapply(1:length(v), function(x){
    getBdryStage2Nparam(ej = gOpt*v[x], aj1= a1[x], ss1 = SS1[x], ss2 = SS2[x])
  }))

  list('gamma'=gOpt, 'PCER_adj'=PCER_adj, 'Stage2AdjBdry'=Stage2AdjBdry)
}




