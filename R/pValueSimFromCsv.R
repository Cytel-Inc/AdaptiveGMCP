# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

#' # The file contains function to run simulation test cases from csv file.----
#' ## Author: Ajoy.M
#'
#' # Run Simulation from CSV ------ -
#' #' @param TestCase CSV with test cases
#' simCSV <- function(TestCases) {
#'   simResults <- data.frame()
#'   for (tc in TestCases$TestCaseID) {
#'     cat("Executing Test Case :", tc, "\n")
#'     simResults <- rbind(simResults, sim1TC(tc))
#'   }
#'   simResults
#' }
#'
#' # Run single test case
#' sim1TC <- function(tc) {
#'   TC <- TestCases[TestCases$TestCaseID == tc, ]
#'   SampleSize <- TC$SampleSize
#'   TailType <- TC$TailType
#'   nArms <- TC$nArms
#'   Arms.Name <- csv2vec(TC$Arms.Name, AsNumeric = F)
#'   Arms.Mean <- csv2vec(TC$Arms.Mean)
#'   Arms.std.dev <- csv2vec(TC$Arms.std.dev)
#'   Arms.alloc.ratio <- csv2vec(TC$Arms.alloc.ratio)
#'   Hypothesis <- csv2vec(TC$Hypothesis, AsNumeric = F)
#'   WI <- csv2vec(TC$WI)
#'   G <- vector2matrix(TC$G, byrow = T, nrow = length(WI))
#'   test.type <- TC$test.type
#'   Correlation <- vector2matrix(TC$Correlation, byrow = T, nrow = length(WI))
#'   info_frac <- csv2vec(TC$info_frac)
#'   Eff_bdry <- csv2vec(TC$Eff_bdry)
#'   MultipleWinners <- TC$MultipleWinners
#'   Selection <- TC$Selection
#'   SelectionLook <- TC$SelectionLook
#'   SelectionMethods <- TC$SelectionMethods
#'   SelectionCriterion <- TC$SelectionCriterion
#'   SelectionParmeter <- csv2vec(TC$SelectionParmeter)
#'   UpdateStrategy <- TC$UpdateStrategy
#'   nSimulation <- TC$nSimulation
#'   Seed <- TC$Seed
#'   ImplicitSSR <- TC$ImplicitSSR
#'   SummaryStat <- TC$SummaryStat
#'
#'   out <- adaptGMCP_SIM(
#'     SampleSize = SampleSize, TailType = TailType, nArms = nArms, Arms.Name = Arms.Name, Arms.Mean = Arms.Mean, Arms.std.dev = Arms.std.dev,
#'     Arms.alloc.ratio = Arms.alloc.ratio, Hypothesis = Hypothesis, WI = WI, G = G, Correlation = Correlation, test.type = test.type,
#'     Eff_bdry = Eff_bdry, info_frac = info_frac, MultipleWinners = MultipleWinners,
#'     Selection = Selection, SelectionLook = SelectionLook, SelectionMethods = SelectionMethods, SelectionCriterion = SelectionCriterion, SelectionParmeter = SelectionParmeter,
#'     UpdateStrategy = UpdateStrategy, ImplicitSSR = ImplicitSSR,
#'     nSimulation = nSimulation, Seed = Seed, SummaryStat = SummaryStat
#'   )
#'   simResultsName <- c(
#'     "TestCaseID", "GlobalPower", "ConjunctivePower",
#'     "DisjunctivePower", "FWER", "Seed", "ElapsedTime"
#'   )
#'   sim1Results <- data.frame(matrix(
#'     c(
#'       tc, unlist(out$Overall_Powers["Values"]),
#'       out$Seed, out$elapsedTime
#'     ),
#'     nrow = 1
#'   ))
#'   names(sim1Results) <- simResultsName
#'   sim1Results
#' }
#'
#' # To convert comma separated inputs to vector
#' csv2vec <- function(inp, AsNumeric = T) {
#'   if (is.numeric(inp) || is.na(inp)) {
#'     return(inp)
#'   } else {
#'     if (AsNumeric) {
#'       return(as.numeric(unlist(strsplit(inp, split = ","))))
#'     } else {
#'       return(unlist(strsplit(inp, split = ",")))
#'     }
#'   }
#' }
#'
#' # To convert comma separated inputs to matrix
#' vector2matrix <- function(inp, byrow, nrow) {
#'   if (is.na(inp)) {
#'     return(inp)
#'   } else {
#'     return(matrix(
#'       as.numeric(sub("NA", NA, unlist(strsplit(inp, split = ",")))),
#'       byrow = byrow, nrow = nrow
#'     ))
#'   }
#' }
