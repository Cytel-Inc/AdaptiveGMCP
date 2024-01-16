#### Compute plan boundary for parametric, non-parametric and mixed case
#' @export
planBdryCER <- function(nHypothesis,nEps,nLooks, alpha,info_frac,typeOfDesign,test.type,Sigma,WH,HypoMap,Scale)
{
  Method <- c() ; SubSets <- c(); SignfLevel <-c()

  Stage1Bdry <- Stage2Bdry <- matrix(0, nrow = nrow(WH),ncol = nHypothesis)
  colnames(Stage1Bdry) <- paste('a', 1:nHypothesis,'1', sep='')
  colnames(Stage2Bdry) <- paste('a', 1:nHypothesis,'2', sep='')

  for(i in 1:nrow(WH))
  {
    #print(i)
    J <- as.numeric(WH[i, 1:nHypothesis])
    w <- as.numeric(WH[i, (nHypothesis+1):(2*nHypothesis)])

    get_Sets <- connSets(J=J,w=w,test.type = test.type, HypoMap=HypoMap)
    conn_Sets <- get_Sets$connSets
    conn_Sets_name <- paste("(", paste(
      unlist(lapply(conn_Sets,function(x){paste(x,collapse = ',')}))
      ,collapse = "),("), ")", sep = "")

    SubSets <- c(SubSets, conn_Sets_name)
    Method <-  c(Method, get_Sets$Method)

    siglev <- list()
    for (edx in conn_Sets) {
      if(length(edx)>1)
      {
        ## Dunnett Weighted Parametric ##
        sig_level <- alpha*sum(w[edx])
        gIDX <- unique(HypoMap[edx,]$Groups)
        wJ <- as.numeric(w) #wJ : weights for the intersection J(including weights with 0)

        plan_parm_bdry <- getPlanParmBdry(gIDX=gIDX, hIDX = edx, alpha = sig_level, nLooks = nLooks,
                                          info_frac = info_frac, wJ = wJ,
                                        Sigma=Sigma, typeOfDesign = typeOfDesign, Scale = Scale)
        Stage1Bdry[i,edx] <- plan_parm_bdry$Stage1Bdry[edx]
        Stage2Bdry[i,edx] <- plan_parm_bdry$Stage2Bdry[edx]

      }else
      {
        ## Non-Parametric ##
        sig_level <- alpha*sum(w[edx])
        plan_nparm_bdry <- getPlanNonParmBdry(nLooks = nLooks, sig_level = sig_level,
                                    info_frac = info_frac, typeOfDesign = typeOfDesign)
        Stage1Bdry[i,edx] <- plan_nparm_bdry$Stage1Bdry
        Stage2Bdry[i,edx] <- plan_nparm_bdry$Stage2Bdry
      }
      siglev <- append(siglev,sig_level)
    }
    SignfLevel <- c(SignfLevel,
                    paste("(", paste(
                      unlist(lapply(siglev,function(x){paste(x,collapse = ',')}))
                      ,collapse = "),("), ")", sep = ""))

  }
  #### Preparation of output tables ###

  #Intersection Weights#
  InterWeightsTab <- WH

  #Test procedure
  TestProcedureTab <- cbind(WH[,1:(ncol(WH)/2)],
                            SubSets,
                            Method,
                            SignfLevel)
  #Stage-1 Boundary
  Stage1BdryTab <- cbind(WH[,1:(ncol(WH)/2)],Stage1Bdry)

  PlanBdryTab <- list('Intersection_Weights'=InterWeightsTab,
                      'Test_Procedure'=TestProcedureTab,
                      'Stage1_Boundary'=Stage1BdryTab)
  if(nLooks == 2)
  {
    Stage2BdryTab <- cbind(WH[,1:(ncol(WH)/2)],Stage2Bdry)
    PlanBdryTab <- list('Intersection_Weights'=InterWeightsTab,
                        'Test_Procedure'=TestProcedureTab,
                        'Stage1_Boundary'=Stage1BdryTab,
                        'Stage2_Boundary'=Stage2BdryTab)
  }

  list('Stage1Bdry'=Stage1Bdry,
       'Stage2Bdry'=Stage2Bdry,
       'PlanBdryTable'=PlanBdryTab)
}








