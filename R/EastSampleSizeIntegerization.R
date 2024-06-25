#Arm-wise look-wise sample size distribution

SS_Alloc=function(SS='SampleSize',r='AllocationRatio',InfoFrac='InfoFrac')
{
  round_off=function(x)
  {
    round_off2=function(x)
    {y=floor(x);   eps=x-y
    abs_diff=abs(eps-0.5)
    if(abs_diff<=1E-12) ## tolarance for absolute diffrences
    {
      eps=0.5
    }else
    {
      eps=eps
    }
    if(eps>=0.5)
    {z=y+1}else{ z=y}
    z }
    sapply(x, round_off2)
  }
  ####################
  narms=length(r)
  nlooks=length(InfoFrac)

  ####################
  f = r/sum(r)                # Allocation Fraction incremental
  f_cum = cumsum(f)           # Allocation Fraction Cumulative
  f_cum[length(f_cum)] = 1

  #########################
  SS_dbl=SS*InfoFrac          # Per look total Sample Size in dbl scale
  SS_int=round_off(SS_dbl)    # Per look total Sample Size in int scale
  SS_int[length(SS_int)]=SS   # Final look SS is the total input sample size
  ##########################

  SS_Arm_Cum=c()              # Arm-wise Sample Size allocation - Cumulative Scale #
  SS_Arm_Inc=c()              # Arm-wise Sample Size allocation - Incremental Scale #

  for(j in 1:nlooks)
  {
    ss_arm_j=round_off(SS_int[j]*f_cum)
    SS_Arm_Cum=rbind(SS_Arm_Cum,ss_arm_j)
  }
  SS_Arm_Inc=cbind(SS_Arm_Cum[,1],t(diff(t(SS_Arm_Cum))))

  df_int=data.frame(cbind(SS_int,SS_Arm_Inc))
  df_int$InfoFrac=df_int$SS_int/SS

  ## the first column is corresponding to look wise total alloted samples and the later columns are representing
  ## the arm-wise(starting from control) distribution of samples for that particular look, the last column is the
  ## recomputed Information fraction
  df_int
}





# Planned Arm-wise & Look-wise samples------------------------- -
getPlanAllocatedSamples <- function(SS, allocRatio, info_frac) {
  col_names <- c("Control", paste("Treatment", 1:(length(allocRatio) - 1), sep = ""))
  SS_Cum_Int_Df <- SS_Alloc(SS = SS, r = allocRatio, InfoFrac = info_frac)
  ss_plan_cum <- SS_Cum_Int_Df[,c(-1,-ncol(SS_Cum_Int_Df))]
  rownames(ss_plan_cum) <- NULL
  colnames(ss_plan_cum) <- col_names

  ss_plan_incr <- ss_plan_cum
  if(length(info_frac) > 1){
    for(lkIDX in 2:length(info_frac)){
      ss_plan_incr[lkIDX, ] <- ss_plan_cum[lkIDX, ]-ss_plan_cum[lkIDX-1, ]
    }
  }
  colnames(ss_plan_cum) <- col_names
  list("IncrementalSamples" = ss_plan_incr, "CumulativeSamples" = ss_plan_cum)
}











