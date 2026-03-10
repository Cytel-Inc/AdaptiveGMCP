# nolint start
#' Generate complete survival trial data upfront for all looks.
#'
#' Generates subject arrival times via piecewise constant accrual,
#' allocates subjects to arms, and generates exponential survival times.
#' Event indicators are NOT computed here; they are determined during
#' look-wise analysis by \code{AnalyzeSurvivalAtLook()}.
#'
#' @param nSubjects Integer. Total sample size N.
#' @param nArms Integer. Number of trial arms (including control).
#' @param Arms.alloc.ratio Numeric vector. Arm-wise allocation ratios
#'   (length = \code{nArms}). First element is control.
#' @param armsHazardRates Numeric vector. Arm-wise hazard rates
#'   (length = \code{nArms}). First element is control arm hazard rate.
#'   All values must be positive.
#' @param accrualStartTimes Numeric vector. Starting time points of
#'   piecewise constant accrual periods. First element must be 0.
#' @param accrualRates Numeric vector. Accrual rates (subjects per unit
#'   time) for each piece. Same length as \code{accrualStartTimes}.
#'   All values must be positive.
#' @param seed Integer or NULL. Random seed for reproducibility.
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{subjectId}{Integer. Subject identifier (1 to N).}
#'     \item{arm}{Integer. Allocated arm (1 = control, 2+ = treatments).}
#'     \item{arrivalTime}{Numeric. Calendar time of subject arrival.}
#'     \item{survivalTime}{Numeric. Time from arrival to event (follow-up scale).}
#'     \item{eventTime}{Numeric. Calendar time of event (arrivalTime + survivalTime).}
#'   }
#'
#' @details
#' Arrival times are generated using piecewise constant accrual with uniform
#' sampling within each piece. For pieces 1 to (n-1), the expected number of
#' subjects is \code{round(accrualRates[i] * duration)}. For the last piece
#' (which has no defined end time), the remaining balance of subjects is
#' enrolled with end time derived as \code{balance / accrualRates[n]}.
#'
#' Treatment allocation is done without replacement so that actual arm-wise
#' counts match the target allocation ratios as closely as possible.
#'
#' Survival times follow an Exponential distribution with arm-specific
#' hazard rates.
#'
#' @keywords internal
GenSurvivalData <- function(nSubjects,
                            nArms,
                            Arms.alloc.ratio,
                            armsHazardRates,
                            accrualStartTimes,
                            accrualRates,
                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # --- 1. Generate arrival times (piecewise constant accrual, uniform) ---
  arrivalTimes <- GenerateArrivalTimes(
    nSubjects = nSubjects,
    accrualStartTimes = accrualStartTimes,
    accrualRates = accrualRates
  )

  # --- 2. Allocate treatments (without replacement) ---
  armAssignments <- AllocateTreatments(
    nSubjects = nSubjects,
    nArms = nArms,
    Arms.alloc.ratio = Arms.alloc.ratio
  )

  # --- 3. Generate survival times (exponential, arm-specific hazard) ---
  survivalTimes <- GenerateSurvivalTimes(
    armAssignments = armAssignments,
    armsHazardRates = armsHazardRates
  )

  # --- 4. Assemble and return data.frame ---
  data.frame(
    subjectId = seq_len(nSubjects),
    arm = armAssignments,
    arrivalTime = arrivalTimes,
    survivalTime = survivalTimes,
    eventTime = arrivalTimes + survivalTimes
  )
}


#' Generate subject arrival times using piecewise constant accrual.
#'
#' For each accrual piece (except the last), subjects are allocated based on
#' rate x duration and sampled uniformly within the piece interval. The last
#' piece (open-ended) enrolls remaining subjects with its end time derived
#' from balance / rate.
#'
#' @param nSubjects Integer. Total number of subjects to enroll.
#' @param accrualStartTimes Numeric vector. Start times of accrual pieces.
#' @param accrualRates Numeric vector. Accrual rates per piece.
#'
#' @return Numeric vector of length \code{nSubjects} with sorted arrival times.
#'
#' @keywords internal
GenerateArrivalTimes <- function(nSubjects,
                                 accrualStartTimes,
                                 accrualRates) {
  nPieces <- length(accrualStartTimes)
  arrivalTimes <- numeric(0)
  subjectsRemaining <- nSubjects

  if (nPieces > 1) {
    for (i in seq_len(nPieces - 1)) {
      if (subjectsRemaining <= 0) break

      pieceDuration <- accrualStartTimes[i + 1] - accrualStartTimes[i]
      expectedSubjects <- round(accrualRates[i] * pieceDuration)
      # Do not exceed remaining subjects
      nInPiece <- min(expectedSubjects, subjectsRemaining)

      if (nInPiece > 0) {
        pieceTimes <- stats::runif(
          n = nInPiece,
          min = accrualStartTimes[i],
          max = accrualStartTimes[i + 1]
        )
        arrivalTimes <- c(arrivalTimes, pieceTimes)
        subjectsRemaining <- subjectsRemaining - nInPiece
      }
    }
  }

  # Last piece: open-ended, derive end time from balance / rate
  if (subjectsRemaining > 0) {
    lastPieceIdx <- nPieces
    lastPieceDuration <- subjectsRemaining / accrualRates[lastPieceIdx]
    lastPieceEnd <- accrualStartTimes[lastPieceIdx] + lastPieceDuration

    lastPieceTimes <- stats::runif(
      n = subjectsRemaining,
      min = accrualStartTimes[lastPieceIdx],
      max = lastPieceEnd
    )
    arrivalTimes <- c(arrivalTimes, lastPieceTimes)
  }

  # Sort chronologically
  sort(arrivalTimes)
}


#' Allocate subjects to treatment arms without replacement.
#'
#' Computes arm-wise sample sizes from allocation ratios, ensuring the total
#' equals \code{nSubjects}, then randomly permutes the assignments.
#'
#' @param nSubjects Integer. Total number of subjects.
#' @param nArms Integer. Number of arms.
#' @param Arms.alloc.ratio Numeric vector. Allocation ratios (length = nArms).
#'
#' @return Integer vector of length \code{nSubjects} with arm assignments
#'   (1 = control, 2+ = treatments).
#'
#' @keywords internal
AllocateTreatments <- function(nSubjects,
                               nArms,
                               Arms.alloc.ratio) {
  # Normalize ratios
  ratioSum <- sum(Arms.alloc.ratio)
  nPerArm <- round(nSubjects * Arms.alloc.ratio / ratioSum)

  # Adjust to ensure total equals nSubjects
  diff <- nSubjects - sum(nPerArm)
  if (diff != 0) {
    # Distribute the difference to the arm(s) with the largest fractional part
    fractional <- (nSubjects * Arms.alloc.ratio / ratioSum) - floor(nSubjects * Arms.alloc.ratio / ratioSum)
    if (diff > 0) {
      # Need to add subjects
      addIdx <- order(fractional, decreasing = TRUE)[seq_len(abs(diff))]
      nPerArm[addIdx] <- nPerArm[addIdx] + 1
    } else {
      # Need to remove subjects
      removeIdx <- order(fractional, decreasing = FALSE)[seq_len(abs(diff))]
      nPerArm[removeIdx] <- nPerArm[removeIdx] - 1
    }
  }

  # Create arm assignment vector and randomly permute
  armAssignments <- rep(seq_len(nArms), times = nPerArm)
  sample(armAssignments)
}


#' Generate exponential survival times based on arm-specific hazard rates.
#'
#' @param armAssignments Integer vector. Arm assignment per subject.
#' @param armsHazardRates Numeric vector. Hazard rate per arm.
#'
#' @return Numeric vector of survival times (follow-up time scale).
#'
#' @keywords internal
GenerateSurvivalTimes <- function(armAssignments,
                                  armsHazardRates) {
  stats::rexp(
    n = length(armAssignments),
    rate = armsHazardRates[armAssignments]
  )
}


#' Determine calendar times at which each look occurs.
#'
#' Uses event-based information fractions to find the calendar time at which
#' the required number of events has accumulated for each look.
#'
#' @param survivalData Data.frame returned by \code{GenSurvivalData()}.
#' @param totalEvents Integer. Total number of events (D) for the trial.
#' @param info_frac Numeric vector. Cumulative information fractions
#'   (length = number of looks). Last element should be 1.
#'
#' @return Numeric vector of calendar times (one per look).
#'
#' @keywords internal
DetermineLookTiming <- function(survivalData,
                                totalEvents,
                                info_frac) {
  # Events per look (cumulative)
  events.per.look <- round(totalEvents * info_frac)
  # Ensure last look uses exactly totalEvents
  events.per.look[length(events.per.look)] <- totalEvents

  # Sort all event times chronologically
  sorted.event.times <- sort(survivalData$eventTime)

  # Calendar time for each look = time of the k-th event
  # TODO: The following code only handles the success path where there are enough events
  # in the data to meet the required number of events per look.
  # We need to think of and handle all the edge cases too.
  vapply(events.per.look, function(k) {
    sorted.event.times[k]
  }, numeric(1))
}


#' Analyze survival data at a specific look using the log-rank test.
#'
#' Applies administrative censoring at the look calendar time, then computes
#' a stratified (by available arms) log-rank Z-statistic and one-sided p-value
#' for each treatment-vs-control hypothesis. Returns a \code{data.frame}
#' matching the structure of \code{getPerLookTestStat()}.
#'
#' @param simID Integer. Current simulation ID.
#' @param lookID Integer. Current look number.
#' @param survivalData Data.frame returned by \code{GenSurvivalData()}.
#' @param lookCalendarTime Numeric. Calendar time for this look.
#' @param ArmsPresent Logical vector (length = nArms). Which arms are still
#'   active.
#' @param HypoMap Data.frame with columns Hypothesis, Groups, EpType,
#'   Control, Treatment.
#'
#' @return A 1-row \code{data.frame} with columns: SimID, LookID,
#'   Delta1..N, StdError1..N, TestStat1..N, RawPvalues1..N, where
#'   N = \code{nrow(HypoMap)}.
#'
#' @details
#' For each hypothesis \eqn{H_i} comparing treatment arm \eqn{k} vs control:
#' \itemize{
#'   \item Subset data to control and treatment \eqn{k} subjects who arrived
#'     before \code{lookCalendarTime}.
#'   \item Apply administrative censoring: event if
#'     \code{eventTime <= lookCalendarTime}, else censored at
#'     \code{lookCalendarTime - arrivalTime}.
#'   \item Compute the log-rank test using \code{survival::survdiff()}.
#'   \item Extract the signed Z-statistic (positive = treatment better than
#'     control, i.e., lower hazard).
#'   \item One-sided p-value: \code{1 - pnorm(Z)}.
#'   \item Delta = log(HR) with SE from the log-rank variance.
#' }
#' Hypotheses for dropped arms are set to \code{NA}.
#'
#' @keywords internal
AnalyzeSurvivalAtLook <- function(simID,
                                  lookID,
                                  survivalData,
                                  lookCalendarTime,
                                  ArmsPresent,
                                  HypoMap) {
  n.hypo <- nrow(HypoMap)

  delta <- se <- test.stat <- p.val <- rep(NA_real_, n.hypo)

  # Create analysis dataset: only subjects who arrived before the look
  analysis.data <- survivalData[survivalData$arrivalTime <= lookCalendarTime, ]
  # Event indicator and observed time
  analysis.data$event <- as.integer(analysis.data$eventTime <= lookCalendarTime)
  analysis.data$obs.time <- ifelse(
    analysis.data$event == 1,
    analysis.data$survivalTime,
    lookCalendarTime - analysis.data$arrivalTime
  )

  for (h.idx in seq_len(n.hypo)) {
    ctr.idx <- HypoMap$Control[h.idx]
    trt.idx <- HypoMap$Treatment[h.idx]

    # Skip if either arm has been dropped
    if (!ArmsPresent[ctr.idx] || !ArmsPresent[trt.idx]) next

    # Subset to control and treatment arms
    pair.data <- analysis.data[analysis.data$arm %in% c(ctr.idx, trt.idx), ]
    # Create group indicator: 0 = control, 1 = treatment
    pair.data$group <- as.integer(pair.data$arm == trt.idx)

    # Need at least 1 event in each group to compute log-rank
    events.per.group <- tapply(pair.data$event, pair.data$group, sum)
    if (any(is.na(events.per.group)) || any(events.per.group == 0)) {
      # Cannot compute test stat — set conservative values
      delta[h.idx] <- 0
      se[h.idx] <- NA_real_
      test.stat[h.idx] <- 0
      p.val[h.idx] <- 1
      next
    }

    # Log-rank test via survival::survdiff
    surv.obj <- survival::Surv(time = pair.data$obs.time, event = pair.data$event)
    lr.fit <- survival::survdiff(surv.obj ~ group, data = pair.data)

    # Extract O-E and variance for the treatment group (group index 2)
    obs.trt <- lr.fit$obs[2]
    exp.trt <- lr.fit$exp[2]
    var.trt <- lr.fit$var[2, 2]

    # Signed Z: positive when treatment has fewer events than expected
    # (i.e., treatment is better = lower hazard)
    z.stat <- (exp.trt - obs.trt) / sqrt(var.trt)

    # Log hazard ratio estimate and SE from log-rank
    log.hr <- (obs.trt - exp.trt) / var.trt
    se.log.hr <- 1 / sqrt(var.trt)

    delta[h.idx] <- log.hr
    se[h.idx] <- se.log.hr
    test.stat[h.idx] <- z.stat
    p.val[h.idx] <- 1 - stats::pnorm(z.stat)
  }

  # Build output matching getPerLookTestStat() format
  sumstatdf <- data.frame(matrix(
    c(simID, lookID, delta, se, test.stat, p.val),
    nrow = 1
  ), row.names = NULL)
  colnames(sumstatdf) <- c(
    "SimID", "LookID",
    paste0("Delta", seq_len(n.hypo)),
    paste0("StdError", seq_len(n.hypo)),
    paste0("TestStat", seq_len(n.hypo)),
    paste0("RawPvalues", seq_len(n.hypo))
  )

  sumstatdf
}
# nolint end