# --------------------------------------------------------------------------------------------------
#
# ©2025 Cytel, Inc.  All rights reserved.  Licensed pursuant to the GNU General Public License v3.0.
#
# --------------------------------------------------------------------------------------------------

library(rhandsontable)

# Module UI for rhandsontable
rhandsontableUI <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  ui <- fluidRow(
    # actionButton(ns("addRow"), "Add Endpoint"),
    rHandsontableOutput(ns("tabularInput")),
    # actionButton(ns("submit"), "Submit", class = "btn-primary")
  )
  return(ui)
}

rhandsontableServer <- function(id, initialData, digits = 4) {
  moduleServer(id, function(input, output, session) {
    values <- reactiveValues(data = NULL)

    observe({
      values$data <- initialData()
    })

    output$tabularInput <- renderRHandsontable({
      if (!is.null(values$data)) {
        rhandsontable(values$data, digits = digits)%>%
          hot_table(rowHeaderWidth = 100,stretchH = "all")
      }
    })
  })
}

# function to convert the output into the required format
processData_tabularInput <- function(inputElement) {
  # Convert input from the Shiny UI to R object
  df <- hot_to_r(inputElement)

  # Initialize an empty list in case df is NULL or has no rows
  myList <- list()

  if (!is.null(df) && nrow(df) > 0) {
    # Generate a list of means
    myList <- lapply(seq_len(nrow(df)), function(i) as.vector(unlist(df[i, ], use.names = FALSE)))
    names(myList) <- rownames(df)
  }

  return(myList)
}
