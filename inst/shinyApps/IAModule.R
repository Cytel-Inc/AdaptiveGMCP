# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

# Define UI for the data input module
IADataInputUI <- function(id) {
  ns <- NS(id)
  rHandsontableOutput(ns("dataInput"))
}

# Define server logic for the data input module
IADataInputServer <- function(id, numRows) {
  moduleServer(id, function(input, output, session) {
    # Reactive expression to initialize or update the dataframe
    dfReactive <- reactive({
      if(is.null(input$dataInput) || nrow(hot_to_r(input$dataInput)) != numRows()) {
        # Initialize or re-initialize the dataframe with the specified number of rows
        initialLook <- seq_len(numRows())
        # Calculate equally spaced InfoFraction values
        infoFractionValues <- seq(from = 1/numRows(), to = 1, length.out = numRows())
        df <- data.frame(Look = initialLook, InfoFraction = infoFractionValues)
      } else {
        df <- hot_to_r(input$dataInput)
        # If numRows has changed, adjust the dataframe to match the new number of rows
        if (nrow(df) != numRows()) {
          initialLook <- seq_len(numRows())
          # Recalculate equally spaced InfoFraction values
          infoFractionValues <- seq(from = 1/numRows(), to = 1, length.out = numRows())
          df <- data.frame(Look = initialLook, InfoFraction = infoFractionValues)
        } else {
          # Update the 'Look' column to be a sequence from 1 to the number of rows
          df$Look <- seq_len(nrow(df))
          # Recalculate InfoFraction to ensure they are equally spaced
          df$InfoFraction <- seq(from = 1/nrow(df), to = 1, length.out = nrow(df))
        }
      }
      return(df)
    })

    # Use the reactive dataframe for rendering the table
    output$dataInput <- renderRHandsontable({
      rh <- rhandsontable(dfReactive(), rowHeaders = FALSE)
      rh <- hot_col(rh, "Look", readOnly = TRUE) # Make the 'Look' column read-only
      rh <- hot_col(rh, "InfoFraction", readOnly = FALSE) # Ensure InfoFraction can be edited if needed
      return(rh)
    })
  })
}

# function to convert the output into the required format
processData_IAInput <- function(inputElement) {
  # Convert input from the Shiny UI to R object
  df <- hot_to_r(inputElement)

  # Initialize an empty list in case df is NULL or has no rows
  myList <- list()

  if (!is.null(df) && nrow(df) > 0) {
    # Generate a list of means
    myList <- lapply(seq_len(nrow(df)), function(i) as.vector(unlist(df[i, 2], use.names = FALSE)))
  }

  return(myList)
}
