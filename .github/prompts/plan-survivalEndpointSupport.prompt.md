# Plan: Add Survival Endpoint Support for P-value Combination Method

## Context

This package has been developed to perform simulations and analysis for clinical trials involving multiple treatment arms and multiple continuous endpoints, or multiple binary endpoints, or multiple mixed (some continuous and some binary) endpoints. It works for both fixed sample (single look) trials as well as trials with 2 looks (one interim look). It implements both p-value combination and CER methods.

The new assignment is to extend this package to simulate trials involving 2 looks (one interim look), a single survival (time to event) endpoint (not multiple endpoints), and multiple treatment arms using the p-value combination method only, and not the CER method. No modification is to be done to the analysis functions.

**CRITICAL CONSTRAINT:** No modifications are to be done to binary and continuous endpoints related code. That code should continue to work exactly as it works today. All survival-specific logic must be implemented in new functions and as branching logic that only executes for survival trials.

## Architecture Overview

The AdaptiveGMCP package uses a type-based branching architecture where endpoint types control data generation, test statistics, and boundaries. Adding survival endpoints requires extending this branching pattern in 4-5 key functions without modifying the analysis functions (which already handle raw p-values generically).

**Key architectural patterns:**
- Endpoint types (`lEpType`) are defined as list parameter (e.g., `list("EP1" = "Continuous", "EP2" = "Binary")`)
- Type-specific branching in data generation, test statistics, and boundary computation
- Raw p-values flow through existing combination machinery without changes needed

### Critical Architectural Difference for Survival Endpoints

**Current incremental data generation (Continuous/Binary):**

`SingleSimCombPValue()` has a while loop over looks:
1. Calls `getInterimSSIncr()` to calculate arm-wise sample sizes for current look
2. Calls `genIncrLookSummary()` for current look, which:
   - Loops over all arms present at current stage
   - Calls `genNormToOther2()` for each arm to generate subject responses
   - Generates only response data (no arrival times, no dropout times)
   - Arm allocation already done by `getInterimSSIncr()`

**Required approach for Survival endpoints:**

This incremental look-by-look data generation **will not work** for survival endpoints because:
- Survival data must be generated **upfront for the entire trial** (combined data for all looks)
- Subject **arrival times** must be generated first (accrual process)
- Then subject **survival times** are generated
- Then the complete dataset is **analyzed look-wise based on calendar time**

**Implications:**
- Need separate data generation path for survival trials
- Generate full trial dataset at simulation start
- Perform look-wise analyses on the same complete dataset
- Look timing determined by accumulated events, not calendar time or sample size

### Specific Implementation Details for Survival Endpoints

**Look timing definition (event-based):**
- Accept total number of events `D` as input (in addition to sample size `N` via existing `SampleSize` parameter)
- `D <= N` (total events cannot exceed total subjects)
- Use existing `info_frac` parameter to define look-wise events: `D * info_frac` (rounded to nearest integer)
- At Look 1 (interim): Count events chronologically until reaching expected events for Look 1
  - Calendar time of last included event = Look 1 calendar time
  - Remaining subjects without events are administratively censored at this time
- At Look 2 (final): Use all subject data
  - Calendar time of last event in trial = total study duration
- Information fractions are event-based, not calendar-time-based

**Arrival time generation (piecewise constant accrual):**
- Input: Two vectors of equal length (inferred, not specified separately)
  - `accrualStartTimes`: Starting time points of accrual pieces (first element = 0)
  - `accrualRates`: Accrual rate (subjects per unit time) for each piece (must be positive)
- Example: `accrualStartTimes = c(0, 10)`, `accrualRates = c(15, 5)`
  - Piece 1: time 0 to 10, 15 subjects/time unit
  - Piece 2: time 10 onwards, 5 subjects/time unit
- Generate arrival times for `N` subjects (from `SampleSize` parameter) using uniform sampling within each piece:
  - For each piece, calculate expected subjects based on rate and duration
  - Allocate subjects to pieces proportional to their capacity
  - Sample arrival times uniformly within each piece
- Accrual stops when `N` subjects are enrolled
- Arrival time of last subject = total accrual duration
- Arrival times independent of arm allocation

**Treatment allocation:**
- Use existing `Arms.alloc.ratio` parameter
- After generating `N` arrival times, randomly assign each subject to an arm
- Allocation without replacement to match target ratios:
  - Calculate arm-wise sample sizes: `N * Arms.alloc.ratio` (rounded, ensuring sum = N)
  - Randomly allocate subjects to achieve these exact counts

**Survival time generation:**
- Generate on follow-up time scale (time from arrival to event)
- Exponential distribution with arm-specific hazard rates
- Input: Vector of hazard rates (length = `nArms`, first element = control)
  - Element 1: control arm hazard rate
  - Element 2: treatment arm 1 hazard rate
  - Element 3: treatment arm 2 hazard rate, etc.
- All hazard rates must be positive
- Event calendar time = arrival time + survival time

**Event indicators:**
- NOT generated upfront with data
- Calculated during look-wise analysis:
  - Look 1: Event if event calendar time <= Look 1 calendar time, else censored
  - Look 2: All remaining subjects have events by definition (Look 2 time = last event time)

## Implementation Steps

### 1. Add survival-specific parameters

**File:** [R/MAMSMEP_SIMULATION_MAIN.R](R/MAMSMEP_SIMULATION_MAIN.R)

Add new parameters to `simMAMSMEP()` function signature:
- `totalEvents` - Total number of events to observe in trial (D), default = NA. Must satisfy `D <= SampleSize` for survival endpoints
- `armsHazardRates` - Numeric vector of arm-wise hazard rates (length = `nArms`), default = NA
  - Element 1: control arm hazard rate
  - Elements 2+: treatment arms hazard rates
  - All values must be positive
- `accrualStartTimes` - Numeric vector of starting times for accrual pieces (first element must be 0), default = c(0)
- `accrualRates` - Numeric vector of accrual rates for each piece (subjects per unit time), default = c(Inf)
  - Must be same length as `accrualStartTimes`
  - All values must be positive
- Note: `TestStatSurv` not needed initially since only log-rank will be implemented

**File:** [R/inputValidationNew.R](R/inputValidationNew.R#L46-L62)

Add "Survival" to allowed `lEpType` values in `valInpsimMAMSMEP()` and add validation:
- If any endpoint is "Survival", validate:
  - `nEps` must equal 1 (single survival endpoint only)
  - `totalEvents` must be specified and satisfy `totalEvents <= SampleSize`
  - `armsHazardRates` must be specified with length = `nArms`, all positive
  - `accrualStartTimes` and `accrualRates` must be specified with equal lengths
  - `accrualStartTimes[1]` must equal 0
  - All `accrualRates` must be positive
  - `nLooks` must equal 2 (only 2-look designs supported)
  - `Method` must equal "CombPValue" (CER not supported for survival)

### 2. Implement survival-specific data generation

**CRITICAL:** Survival data generation requires a fundamentally different approach than continuous/binary.

**New function needed:** Create `GenSurvivalData()` in new file **R/survivalAnalysis.R**

**Function signature:**
```r
GenSurvivalData <- function(
  nSubjects,           # Total sample size N (from SampleSize parameter)
  nArms,               # Number of arms
  armsAllocRatio,      # Arm-wise allocation ratios
  armsHazardRates,     # Arm-wise hazard rates (length = nArms)
  accrualStartTimes,   # Starting times of accrual pieces
  accrualRates,        # Accrual rates for each piece
  seed                 # Random seed
)
```

**Function logic:**
1. **Generate arrival times** using piecewise constant accrual with uniform sampling:
   - For pieces 1 to (n-1) where n = number of accrual pieces:
     - Calculate duration of each piece: `duration = accrualStartTimes[i+1] - accrualStartTimes[i]`
     - Calculate expected subjects for piece: `round(accrualRates[i] × duration)`
     - Sample arrival times uniformly: `Uniform(accrualStartTimes[i], accrualStartTimes[i+1])`
   - For the last piece (piece n, which has no defined end time):
     - Calculate remaining subjects: `balance = N - (subjects allocated to pieces 1 to n-1)`
     - Calculate expected duration: `duration = balance / accrualRates[n]`
     - Calculate end time: `endTime = accrualStartTimes[n] + duration`
     - Sample arrival times uniformly: `Uniform(accrualStartTimes[n], endTime)`
   - Combine and sort all arrival times chronologically
   - Arrival time of last subject = total accrual duration
   - Result: Vector of N arrival times in calendar time

2. **Allocate treatments** (without replacement):
   - Calculate arm-wise sample sizes: `nPerArm = round(nSubjects * Arms.alloc.ratio / sum(Arms.alloc.ratio))`
   - Adjust to ensure `sum(nPerArm) == nSubjects`
   - Randomly permute arm assignments to N subjects

3. **Generate survival times** (follow-up time scale):
   - For each subject with assigned arm k:
     - Generate survival time ~ Exponential(rate = armsHazardRates[k])
   - Calculate event calendar time = arrival time + survival time

4. **Return data.frame** with columns:
   - `subjectId`: 1 to N
   - `arm`: Allocated arm (1 = control, 2+ = treatments)
   - `arrivalTime`: Calendar time of arrival
   - `survivalTime`: Time from arrival to event (follow-up scale)
   - `eventTime`: Calendar time of event (arrivalTime + survivalTime)
   - **DO NOT include event indicator** - computed during analysis

**Call location:** This function should be called **once at the start** of `SingleSimCombPValue()`, **before** the while loop over looks.

**Do NOT extend `genNormToOther2()`** - this function is designed for incremental look-wise generation and is not suitable for survival endpoints.

### 3. Implement survival-specific test statistics and look determination

**New function needed:** Create `DetermineLookTiming()` in **R/survivalAnalysis.R**

**Function signature:**
```r
DetermineLookTiming <- function(
  survivalData,        # Complete survival dataset
  totalEvents,         # D = total events
  info_frac            # Information fractions (length = 2) - existing parameter
)
```

**Function logic:**
1. Calculate expected events per look: `eventsPerLook = round(totalEvents * info_frac)`, ensuring total doesn't exceed D
2. Sort all event times in `survivalData$eventTime` chronologically
3. Look 1 calendar time = eventTime of the `eventsPerLook[1]`-th event
4. Look 2 calendar time = eventTime of the last event (maximum eventTime)
5. Return vector of look calendar times

**New function needed:** Create `AnalyzeSurvivalAtLook()` in **R/survivalAnalysis.R**

**Function signature:**
```r
AnalyzeSurvivalAtLook <- function(
  survivalData,        # Complete survival dataset
  lookCalendarTime,    # Calendar time for this look
  ArmsPresent,         # Vector of arms still in trial - existing variable
  nArms,               # Total number of arms
  nHypothesis          # Total number of hypotheses
)
```

**Function logic:**
1. **Create analysis dataset** for current look:
   - For each subject:
     - If `eventTime <= lookCalendarTime`: Include as event (status = 1, time = survivalTime)
     - If `eventTime > lookCalendarTime`: Administratively censor (status = 0, time = lookCalendarTime - arrivalTime)
   - Filter to include only subjects with `arrivalTime <= lookCalendarTime`

2. **Compute log-rank test statistics** for each hypothesis (treatment vs control):
   - For each treatment arm k in ArmsPresent (k >= 2):
     - Subset data to control (arm=1) and treatment k (arm=k)
     - Compute log-rank test using survival package or manual calculation
     - Extract Z-statistic (signed for direction), p-value (one-sided)
     - Estimate log hazard ratio (delta) and its SE if needed

3. **Return structure** matching what `getPerLookTestStat()` returns:
   - Named list with elements for each hypothesis
   - Each element contains: `delta`, `se`, `teststat`, `pvalue`
   - Set NA for hypotheses of dropped arms

**Do NOT modify `getPerLookTestStat()`** in [R/genBinRespMAMSMEP.R](R/genBinRespMAMSMEP.R) since this is designed for incremental look data, not the upfront survival data structure.

**Integration point:** In `SingleSimCombPValue()`, add branching logic:
```r
if (isSurvivalTrial) {
  if (currentLook == 1) {
    lookCalendarTimes <- determineLookTiming(survivalData, TotalEvents, info_frac)
  }
  testResults <- analyzeSurvivalAtLook(
    survivalData, lookCalendarTimes[currentLook], ArmsPresent, nArms, nHypothesis
  )
} else {
  testResults <- getPerLookTestStat(incrementalData, ...)
}
```

### 4. Update boundary computation for survival

**File:** [R/CER_boundary.R](R/CER_boundary.R#L240-L320)

**Note:** Since CER method is NOT being implemented for survival (only p-value combination), boundary computation changes are **NOT NEEDED**.

For p-value combination method, boundaries come from `rpact::getDesignGroupSequential()` which provides alpha spending boundaries based on `info_frac` parameter. Since information fractions for survival are event-based (not calendar time), the existing `info_frac` parameter directly represents the proportion of total events at each look.

No modifications to `getSigma()` or boundary computation functions are required for survival endpoints with p-value combination method.

### 5. Modify simulation flow for survival trials

**File:** Need to identify where `SingleSimCombPValue()` is defined (likely in R/MAMSMEP_SIMULATION_MAIN.R or related file)

**Key changes:**

1. **Add survival detection logic** at the start of the function:
```r
isSurvivalTrial <- any(sapply(lEpType, function(x) x == "Survival"))
```

2. **Add upfront data generation** before the while loop:
```r
if (isSurvivalTrial) {
  # Generate complete survival dataset for entire trial
  survivalData <- GenSurvivalData(
    nSubjects = SampleSize,
    nArms = nArms,
    armsAllocRatio = Arms.alloc.ratio,      # Existing parameter name
    armsHazardRates = armsHazardRates,      # New parameter name
    accrualStartTimes = accrualStartTimes,  # New parameter name
    accrualRates = accrualRates,            # New parameter name
    seed = currentSeed
  )
  
  # Determine look calendar times based on events
  lookCalendarTimes <- NULL  # Will be computed at first look
}
```

3. **Branch test statistic computation** inside the while loop:
```r
if (isSurvivalTrial) {
  # Determine look timing at first look
  if (currentLook == 1) {
    lookCalendarTimes <- DetermineLookTiming(
      survivalData = survivalData,
      totalEvents = totalEvents,
      info_frac = info_frac               # Existing parameter name
    )
  }
  
  # Analyze survival data at current look
  testResults <- AnalyzeSurvivalAtLook(
    survivalData = survivalData,
    lookCalendarTime = lookCalendarTimes[currentLook],
    ArmsPresent = ArmsPresent,            # Existing variable name
    nArms = nArms,
    nHypothesis = nHypothesis
  )
} else {
  )
} else {
  # Existing incremental approach for continuous/binary
  # Get sample sizes for current look
  lookSSInfo <- getInterimSSIncr(...)
  
  # Generate incremental data
  lookSummary <- genIncrLookSummary(...)
  
  # Compute test statistics
  testResults <- getPerLookTestStat(...)
}
```

4. **Note on sample size functions:** For survival trials, `getInterimSSIncr()` and `genIncrLookSummary()` are NOT called. All data is pre-generated and look-wise analysis uses the same complete dataset.

### 6. Update helper functions and UI

**File:** [inst/shinyApps/AdaptGMCPSimApp.R](inst/shinyApps/AdaptGMCPSimApp.R)

Extend UI modules to accept survival parameters when endpoint type is "Survival":
- Input fields for hazard ratios
- Input fields for control event rate
- Input fields for accrual duration
- Input fields for follow-up duration
- Test statistic selection (LogRank/Cox)

**Note:** `armSumry()` in [R/genBinRespMAMSMEP.R](R/genBinRespMAMSMEP.R) may not need modification if it's only used for incremental data. Survival summaries (events, person-time) will be computed within `analyzeSurvivalAtLook()`.

## Design Decisions to Make

### 1. Information fraction definition for survival
**DECISION MADE:** Event-based information fractions
- `info_frac` represents proportion of total events (D) at each look
- Look timing determined by counting events chronologically
- Existing `info_frac` parameter used without modification

### 2. Look timing specification
**DECISION MADE:** Event-based automatic determination
- Look calendar times are NOT specified by user
- Derived automatically by counting D * info_frac events chronologically
- Look 1 time = calendar time of last event included in Look 1
- Look 2 time = calendar time of last event in entire trial

### 3. Sample size and events handling
**DECISION MADE:** Dual specification
- `SampleSize` parameter = N = total subjects to enroll (existing parameter)
- `TotalEvents` parameter = D = total events expected (new parameter)
- Must satisfy D <= N
- Both specified by user

### 4. Correlation handling
**DECISION MADE:** Ignore and validate
- `EP.Corr` parameter ignored for survival trials
- Validation enforces `nEps = 1` when endpoint type is "Survival"
- Single survival endpoint only in this implementation

### 5. Test statistic choice
**DECISION MADE:** Log-rank only
- Implement log-rank test statistic only
- Cox regression deferred to future enhancement
- No `TestStatSurv` parameter needed initially

### 6. Censoring model
**DECISION MADE:** Administrative censoring only
- Subjects censored at look calendar time if event hasn't occurred yet
- No random censoring (drop-outs) in initial implementation
- Random censoring deferred to future enhancement

### 7. Validation strategy
**DECISION MADE:** Separate test file
- Create `tests/testthat/test-survival.R` for survival-specific tests
- Keep existing continuous/binary tests unchanged

### 8. Treatment selection for survival trials
**DECISION MADE:** P-value based selection only
- Implement treatment selection at interim based on p-values only
- Use existing selection framework (SelectionScale = "pvalue")
- Hazard ratio-based selection deferred to future enhancement

## Key Files to Modify

1. **R/MAMSMEP_SIMULATION_MAIN.R** - Add survival parameters to `simMAMSMEP()`
2. **R/inputValidationNew.R** - Add "Survival" validation, enforce single endpoint for survival
3. **R/survivalAnalysis.R** - **NEW FILE** - Create `genSurvivalData()` and `analyzeSurvivalAtLook()`
4. **R/MAMSMEP_sim_main.R** or equivalent - Modify `SingleSimCombPValue()` to handle survival trials with upfront data generation
5. **R/CER_boundary.R** - Possibly extend `getSigma()` if event-based information fractions are needed
6. **R/simMAMSMEP_Wrapper.R** - Update wrapper to handle survival parameters
7. **inst/shinyApps/AdaptGMCPSimApp.R** - Add UI for survival inputs
8. **tests/testthat/test-survival.R** - **NEW FILE** - Test file for survival endpoints

## Testing Strategy

1. Create simple test case: 3 arms (1 control + 2 treatments), 1 survival endpoint, 2 looks
2. Validate against known results from survival power calculations
3. Test parameter edge cases (high censoring, varying hazard ratios)
4. Verify p-value combination method works correctly with survival p-values
5. Compare with manual log-rank test calculations

## Implementation Order

1. **Phase 1 - Core functionality (P-value combination only)**
   - Add survival parameters to `simMAMSMEP()`:
     - `TotalEvents`, `Arms.HazardRates`, `AccrualStartTimes`, `AccrualRates`
   - Add validation in `valInpsimMAMSMEP()`:
     - Enforce `nEps = 1`, `nLooks = 2`, `Method = "CombPValue"` for survival
     - Validate `TotalEvents <= SampleSize`
     - Validate hazard rates and accrual parameters
   - Create **R/survivalAnalysis.R** with three functions:
     - `GenSurvivalData()`: Piecewise accrual, treatment allocation, exponential survival times
     - `DetermineLookTiming()`: Event-based look timing calculation
     - `AnalyzeSurvivalAtLook()`: Administrative censoring, log-rank test, p-values
   - Modify `SingleSimCombPValue()` to branch for survival trials:
     - Detect survival endpoint via `lEpType`
     - Generate complete data upfront before while loop
     - Skip `getInterimSSIncr()` and `genIncrLookSummary()` for survival
     - Call `DetermineLookTiming()` at first look
     - Call `AnalyzeSurvivalAtLook()` in each look iteration
   - Basic tests with 3-arm, 1-endpoint, 2-look scenario:
     - Verify event counts match expected
     - Verify look timing based on events
     - Verify log-rank p-values are reasonable

2. **Phase 2 - Integration and validation**
   - Update `simMAMSMEP_Wrapper()` to handle survival parameters from CSV
   - Validate against known survival trial power calculations (e.g., rpact)
   - Test edge cases:
     - High censoring scenarios
     - Varying hazard ratios
     - Different allocation ratios
   - Ensure existing continuous/binary tests still pass

3. **Phase 3 - Enhanced features**
   - Implement treatment selection for survival trials (based on p-values)
   - Add Shiny UI components for survival inputs
   - Add Cox regression as alternative to log-rank
   - Add random censoring (exponential dropout)
   - Event-based information fractions (if needed)
   - Comprehensive documentation and vignettes
   - Performance optimization for large sample sizes

4. **Phase 4 - Polish and documentation**
   - Create example scripts in `internalData/`
   - Update package vignettes
   - Add survival-specific validation examples
   - Code review and style compliance (Google R Style Guide)
   - Final integration testing
