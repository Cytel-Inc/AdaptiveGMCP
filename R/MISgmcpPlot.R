# # Function to plot the graph by using initial weights and the transition matrix.
# gmcpPlot <- function(
#     WI = c(1 / 4, 1 / 4, 1 / 4, 1 / 4),
#     G = matrix(
#       c(
#         0, 0.75, 0.25, 0,
#         0.75, 0, 0, 0.25,
#         1, 0, 0, 0,
#         0, 1, 0, 0
#       ),
#       nrow = 4, byrow = TRUE
#     ),
#     nameHypotheses,
#     hGroup,
#     cordinates) {
#   nHypotheses <- length(WI)
#   m <- G
#   ################ Settings ##########################
#   wchar <- "Weight" # Text to appear before weights
#   digits <- 3 # digits to display
#   size <- 4 # Size of the fonts inside nodes
#   boxtextsize <- 3 # Size of the fonts inside transition weights box
#   trhw <- 0.15 # width of the transition weights box
#   trhh <- 0.1 # width of the transition weights box
#   OffSet <- pi / 10 # Distance between two directed lines
#   ###################################################
#
#   # Classes to color code
#   if (!missing(hGroup)) {
#     fill <- as.integer(as.factor(hGroup))
#   } else {
#     fill <- 1:nHypotheses
#   }
#
#   # Color codes
#   color <- generate_colors(length(unique(fill)))
#
#   # manually specifying the node position
#   if (!missing(cordinates)) {
#     x <- sapply(cordinates, "[[", 1)
#     y <- sapply(cordinates, "[[", 2)
#     fig <- gMCPLite::hGraph(
#       nHypotheses = nHypotheses, alphaHypotheses = WI, m = m, wchar = wchar,
#       digits = digits, size = size, boxtextsize = boxtextsize, trhw = trhw, trhh = trhh, offset = OffSet,
#       x = x, y = y,
#       fill = fill,
#       palette = color
#     )
#   } else {
#     fig <- gMCPLite::hGraph(
#       nHypotheses = nHypotheses, alphaHypotheses = WI, m = m, wchar = wchar,
#       digits = digits, size = size, boxtextsize = boxtextsize, trhw = trhw, trhh = trhh, offset = OffSet,
#       fill = fill,
#       palette = color
#     )
#   }
#   return(fig)
# }
