library(shiny)
library(rhandsontable)


transitionMatrixInputUI <- function(id) {
  ns <- NS(id)
  rHandsontableOutput(ns("transitionMatrix"))
}

transitionMatrixInput <- function(input, output, session, dimension) {
  observe({
    n <- as.numeric(dimension())
    
    # Create an n x n matrix
    MAT <- matrix(rep(0, n^2), nrow = n, ncol = n)
    diag(MAT) <- 0
    
    # Set row and column names
    row_names <- paste("H", 1:n) # Creates row names H, H2, ..., Hn
    col_names <- paste("H", 1:n) # Creates column names H1, H2, ..., Hn
    rownames(MAT) <- row_names
    colnames(MAT) <- col_names
    
    output$transitionMatrix <- renderRHandsontable({
      rh_table <- rhandsontable(MAT, readOnly = FALSE) 
      
      # Apply conditional readOnly settings
      for (i in 1:nrow(MAT)) {
        for (j in 1:ncol(MAT)) {
          if (i == j) { # Cells below the diagonal
            rh_table <- hot_cell(rh_table, i, j, readOnly = TRUE)
          }
        }
      }
      
      rh_table %>%
        hot_cols(renderer = "
                  function (instance, td, row, col, prop, value, cellProperties) {
                    Handsontable.renderers.TextRenderer.apply(this, arguments);
                    if (row == col) {
                      td.style.background = 'lightgrey';
                    } else if (col > row) {
                      td.style.background = '';
                    } else {
                      td.style.background = '';
                      Handsontable.renderers.TextRenderer.apply(this, arguments, '', 'NA');
                    }
                  }")
    })
  })
}

processData_transitionMatrix <- function(inputElement) {
  # Convert input from the Shiny UI to R object
  df <- hot_to_r(inputElement)
  
  if (!is.null(df) && nrow(df) > 0) {
    # Generate the matrix
    myMat <- matrix(df, nrow(df), nrow(df))
  }
  return(myMat)
}