# File: sim_graphicalMCP.R
# Utility function for simulations using the graphicalMCP package

library(graphicalMCP)
library(purrr)

#' Function to create a generic graph using the graphicalMCP package
#'
#' @param numHypo number of hypotheses to be tested
#' @param mcp MCP procedure to be used for multiplicity adjustment
#'    mcp values: 1 => Bonferroni, 2 => Weighted Bonferroni, 3 => Holm,
#'    4 => Weighted Holm, 5 => Fixed Sequence, 6 => Fallback, 7 => Hochberg,
#'    8 => Hommel, 9 => Sidak, 10 => Dunnett's Single Step,
#'    11 => Weighted Dunnett's Single Step
#' @param nodeWts weights for the graph nodes specified as a numeric vector
#' @param transMatrix transition matrix for the graph specified as a numeric matrix
#' @param plotGraph if TRUE (default), the graph is plotted
#'
#' Note that you can specify either (numHypo, mcp) OR (nodeWts, transMatrix),
#' but not both. If (numHypo, mcp) are specified, the other parameters are ignored
#'
CreateGenericGraph_gM <- function(numHypo=0, mcp=NULL,
                                  nodeWts=NULL, transMatrix=NULL,
                                  plotGraph=TRUE) {
  # At least some arguments must be set
  stopifnot(!every(c(numHypo, mcp, nodeWts, transMatrix), is.null))

  approach1 <- every(c(numHypo, mcp), is.null)
  approach2 <- every(c(nodeWts, transMatrix), is.null)

  if(approach1 && approach2){
    stop("You must specify either (numHypo, mcp) OR (nodeWts, transMatrix).")
  }

  if(!approach1 && !approach2){
    warning("You have specified both (numHypo, mcp) and (nodeWts, transMatrix).
            Only one of them is required. Ignoring (nodeWts, transMatrix).")
  }

  gr <- NULL

  if(approach1){
    # User would like to use a standard MCP
    if (!(mcp %in% 1:11)){
      stop("Incorrect mcp value specified! Please input a number >= 1 and <= 11.")
    }

    if(!(numHypo > 0 && numHypo <= 25)){
      stop("Incorrect numHypo value specified! Please input a number > 0 and <= 25.")
    }


  }
  else{ # i.e. approach2 is being used

  }

  # Plotting the graph conditionally
  if(plotGraph){
    plot(gr)
  }

  return(gr)
}
