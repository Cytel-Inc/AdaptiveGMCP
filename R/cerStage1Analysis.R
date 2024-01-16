
#Perform Stage-1 Test -------------------------- -
#' @export
PerformStage1Test <-function(
    nArms = 3,
    nEps  = 2,
    nLooks = 2,
    nHypothesis = nEps*(nArms-1),
    sigma = list('Group1' = c(1,1,1), 'Group1' = c(1,1,1)),
    allocRatio = c(1,1,1),
    SampleSize = 500,
    alpha = 0.025,
    info_frac = c(0.5,1),
    typeOfDesign = 'asOF',
    des.type = 'MAMSMEP',
    test.type = 'Partly-Parametric',
    Stage1Pvalues = c(0.00045,0.0952,0.0225,0.1104),
    HypoMap,
    WH
)
{
  #Stage-Wise Cumulative Sample Size
  SS_alloc <- getPlanAllocatedSamples(SS = SampleSize, allocRatio = allocRatio, info_frac = info_frac)
  SS_Cum <- SS_alloc$CumulativeSamples

  #Computed covariance matrix
  if(test.type == 'Partly-Parametric' || test.type == 'Parametric')
  {
    Sigma <- getSigma(SS_Cum = SS_Cum, sigma = sigma, allocRatio = allocRatio)
  }else
  {
    Sigma <- NA
  }

  #Planned Boundaries
  plan_Bdry <- planBdryCER(nHypothesis = nHypothesis,nEps=nEps,nLooks=nLooks,
                           alpha=alpha,info_frac=info_frac,typeOfDesign=typeOfDesign,test.type=test.type,
                           Sigma=Sigma,WH=WH,HypoMap=HypoMap,Scale = 'Score')

  #Stage1 Analysis
  Stage1Analysis <- closedTest(WH = WH,
                               boundary = plan_Bdry$Stage1Bdry,
                               pValues = Stage1Pvalues)

  Stage1Tables <- list('PlannedSampleSize(Cum.)'=SS_Cum,
                  'Stage1_Tables'=plan_Bdry$PlanBdryTable,
                 'Test_Intersection_Hypothesis'=Stage1Analysis$IntersectHypoTest,
                 'Rejection_Status'=Stage1Analysis$PrimaryHypoTest)

  Stage1Obj <- list('HypoMap'= HypoMap, 'info_frac'=info_frac,'AllocSampleSize'=SS_Cum,
                     'Sigma'= Sigma, 'WH' = WH, 'plan_Bdry'=plan_Bdry, 'Stage1Analysis'=Stage1Analysis)

  list('Stage1Tables'=Stage1Tables, 'Stage1Obj'=Stage1Obj)
}





#Planned Arm-wise & Look-wise samples------------------------- -
getPlanAllocatedSamples <- function(SS, allocRatio, info_frac)
{
  col_names <- c('Control',paste('Treatment',1:(length(allocRatio)-1), sep=''))
  ss_frac <- allocRatio/sum(allocRatio)
  ss_lk <- round(SS*info_frac)
  ss_lk[length(info_frac)] <- SS
  ss_incr <- c(ss_lk[1], diff(ss_lk))
  ss_plan_incr <- t(sapply(ss_incr, function(x) round(x*ss_frac)))

  if(length(info_frac)>1)ss_plan_incr[,1] <- ss_incr - rowSums(ss_plan_incr[,-1])

  ss_plan_incr <- data.frame(ss_plan_incr)
  colnames(ss_plan_incr) <- col_names

  #ss_plan_cum <- matrix(NA, nrow = nrow(ss_plan_incr), ncol = ncol(ss_plan_incr))
  ss_plan_cum <- data.frame(matrix(NA, nrow = nrow(ss_plan_incr), ncol = ncol(ss_plan_incr)))
  for(i in 1:nrow(ss_plan_incr)){
    if(i == 1){
      ss_plan_cum[i,]=ss_plan_incr[i,]
    }else
    {
      ss_plan_cum[i,]=ss_plan_incr[i,]+ss_plan_cum[i-1,]
    }
  }
  colnames(ss_plan_cum) <- col_names
  list('IncrementalSamples'=ss_plan_incr, 'CumulativeSamples'=ss_plan_cum)

}


