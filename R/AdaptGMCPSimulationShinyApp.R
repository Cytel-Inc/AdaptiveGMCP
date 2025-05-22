# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

#' Shiny App for AdaptGMCP Simulations
#' @export
AdaptGMCPSimApp <- function() {
  pkg.path <- system.file(package = "AdaptGMCP")
  app.path <- paste0(pkg.path,'/shinyApps/AdaptGMCPSimApp.R')

  tryCatch(
    shiny::runApp(app.path, display.mode = "normal"),
    error = function(err){
      print("Not Able to launch the App")
    }
  )
}
