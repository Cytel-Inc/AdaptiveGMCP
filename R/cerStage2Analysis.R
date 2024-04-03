PerformStage2Test <- function(mcpObj, AdaptStage2) {
  Stage1Objs <- mcpObj$Stage1Obj
  if (!AdaptStage2) {
    WH_modified_idx <- as.vector(apply(
      mcpObj$WH[, grep("H", names(mcpObj$WH))], 1,
      function(x) {
        paste(x, collapse = "")
      }
    ))

    WH_old_idx <- as.vector(apply(
      mcpObj$WH_Prev[, grep("H", names(mcpObj$WH_Prev))], 1,
      function(x) {
        paste(x, collapse = "")
      }
    ))

    oldIdx <- unlist(lapply(
      1:nrow(mcpObj$WH),
      function(i) {
        which(WH_old_idx == WH_modified_idx[i])
      }
    ))

    boundary <- Stage1Objs$plan_Bdry$Stage2Bdry[oldIdx, ]

    Stage2Analysis <- closedTest(
      WH = mcpObj$WH,
      boundary = boundary,
      pValues = mcpObj$p_raw,
      Stage1RejStatus = mcpObj$rej_flag_Prev
    )

    intHypTab <- Stage2Analysis$IntersectHypoTest
    HypoTab <- intHypTab[, grep("H", names(intHypTab))]
    InterHyp <- apply(HypoTab, 1, function(h) {
      paste(names(HypoTab)[which(h == 1)], collapse = ",")
    })

    intRejStat <- sapply(intHypTab$Rejected, function(x) {
      ifelse(x, "Rejected", "Not_Rejected")
    })
    IntersectHypoTest <- knitr::kable(data.frame(
      "Hypotheses" = InterHyp,
      "Status" = intRejStat,
      row.names = NULL
    ), align = "c")

    rejTab <- Stage2Analysis$PrimaryHypoTest
    FinRejStat <- sapply(rejTab, function(x) {
      ifelse(x, "Rejected", "Not_Rejected")
    })
    FinalRejTab <- knitr::kable(data.frame(
      "Hypotheses" = names(rejTab),
      "Status" = FinRejStat, row.names = NULL
    ), align = "c")

    Stage2Tables <- list(
      "Test_Intersection_Hypothesis" = IntersectHypoTest,
      "Final_Rejection_Status" = FinalRejTab
    )
  } else {
    Stage2Analysis <- closedTest(
      WH = mcpObj$WH,
      boundary = mcpObj$AdaptObj$Stage2AdjBdry,
      pValues = mcpObj$p_raw,
      Stage1RejStatus = mcpObj$rej_flag_Prev
    )
    intHypTab <- Stage2Analysis$IntersectHypoTest
    HypoTab <- intHypTab[, grep("H", names(intHypTab))]
    InterHyp <- apply(HypoTab, 1, function(h) {
      paste(names(HypoTab)[which(h == 1)], collapse = ",")
    })

    intRejStat <- sapply(intHypTab$Rejected, function(x) {
      ifelse(x, "Rejected", "Not_Rejected")
    })
    IntersectHypoTest <- knitr::kable(data.frame(
      "Hypotheses" = InterHyp,
      "Status" = intRejStat,
      row.names = NULL
    ), align = "c")

    rejTab <- Stage2Analysis$PrimaryHypoTest
    FinRejStat <- sapply(rejTab, function(x) {
      ifelse(x, "Rejected", "Not_Rejected")
    })
    FinalRejTab <- knitr::kable(data.frame(
      "Hypotheses" = names(rejTab),
      "Status" = FinRejStat, row.names = NULL
    ), align = "c")

    Stage2Tables <- list(
      "Adapt_Test_Tables" = mcpObj$AdaptObj$Stage2Tables,
      "Test_Intersection_Hypothesis" = IntersectHypoTest,
      "Final_Rejection_Status" = FinalRejTab
    )
  }

  list("Stage2Tables" = Stage2Tables, "RejStat" = Stage2Analysis$PrimaryHypoTest)
}

# To modify the stage-2 sample size(SSR)
do_ModifyStage2Sample <- function(allocRatio, ArmsPresent, AllocSampleSize) {
  ArmsPresent <- sort(ArmsPresent, decreasing = F)
  newAllocSampleSize <- AllocSampleSize
  newallocRatio <- allocRatio

  newAllocSampleSize[2, which(!(1:ncol(newAllocSampleSize) %in% ArmsPresent))] <-
    newallocRatio[which(!(1:ncol(newAllocSampleSize) %in% ArmsPresent))] <- NA

  SSRFlag <- readline(prompt = paste("Modify Stage-2 cumulative Sample Size (y/n) :"))
  if (SSRFlag == "y") {
    availArmsName <- names(AllocSampleSize)[ArmsPresent]
    cat("Planned Sample Size for reference: \n")
    print(AllocSampleSize)

    cat(
      "Enter Stage-2 sample size for ", paste(availArmsName, collapse = ", "),
      "(e.g. 100,150)", "\n"
    )
    newSS <- readline()

    newAllocSampleSize[2, ArmsPresent] <- as.numeric(stringr::str_trim(
      unlist(strsplit(newSS, split = ",")),
      "both"
    ))
    newallocRatio <- as.numeric(newAllocSampleSize[2, ]) / as.numeric(newAllocSampleSize[2, 1])

    list(
      "newAllocSampleSize" = newAllocSampleSize,
      "newallocRatio" = newallocRatio
    )
  } else {
    list(
      "newAllocSampleSize" = AllocSampleSize,
      "newallocRatio" = allocRatio
    )
  }
}
