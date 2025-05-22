# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

library(shiny)
library(rhandsontable)


matrixInputUI <- function(id) {
  ns <- NS(id)
  rHandsontableOutput(ns("corrMatrix"))
}

matrixInput <- function(input, output, session, dimension, simulateTrigger) {
  # Reactive value to store the matrix
  matrixData <- reactiveVal()

  observe({
    n <- as.numeric(dimension())

    # Create an n x n matrix
    MAT <- matrix(rep(0, n^2), nrow = n, ncol = n)
    diag(MAT) <- 1
    MAT[upper.tri(MAT)] <- 0
    MAT[lower.tri(MAT)] <- t(MAT)[upper.tri(MAT)]

    rownames(MAT) <- paste0("EP", seq_len(nrow(MAT)))
    colnames(MAT) <- paste0("EP", seq_len(nrow(MAT)))

    # Store the initial matrix in reactive value
    matrixData(MAT)

    output$corrMatrix <- renderRHandsontable({
      rh_table <- rhandsontable(matrixData(), readOnly = FALSE)

      # Apply conditional readOnly settings
      for (i in 1:nrow(MAT)) {
        for (j in 1:ncol(MAT)) {
          if (i > j) { # Cells below the diagonal
            rh_table <- hot_cell(rh_table, i, j, readOnly = TRUE)
          }
        }
      }

      rh_table %>%
        hot_cols(renderer = "
                  function (instance, td, row, col, prop, value, cellProperties) {
                    Handsontable.renderers.TextRenderer.apply(this, arguments);
                    if (row == col) {
                      td.style.background = 'skyblue';
                    } else if (col > row) {
                      td.style.background = '';
                    } else {
                      td.style.background = 'lightgrey';
                      Handsontable.renderers.TextRenderer.apply(this, arguments, '', 'NA');
                    }
                  }")
    })
  })
  # Update reactive matrix on change
  observeEvent(input$corrMatrix, {
    newMat <- hot_to_r(input$corrMatrix)
    if (!is.null(newMat)) {
      symMat <- matrixData()
      symMat[upper.tri(symMat)] <- t(newMat)[upper.tri(newMat)]
      matrixData(symMat)
    }
  })
}

processData_corrMatrix <- function(inputElement) {
  # Convert input from the Shiny UI to R object
  df <- hot_to_r(inputElement)

  if (!is.null(df) && nrow(df) > 0) {
    # Generate the symmetric matrix
    myMat <- matrix(df, nrow(df), nrow(df))
    myMat[lower.tri(myMat)] <- myMat[upper.tri(myMat)]
  }
  return(myMat)
}
