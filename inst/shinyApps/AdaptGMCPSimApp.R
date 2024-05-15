# FIXME number of looks only 1 and 2. multiple winners, selection and implicit ssr only available for number of looks = 2

library(shiny)
library(AdaptGMCP)
library(rhandsontable)
library(future)
library(promises)
library(dplyr)

source("tabularModule.R")
source("correlationMatrixModule.R")
source("transitionMatrixModule.R")
source("IAModule.R")
source("helper_inputDataCsv.R")

ui <- fluidPage(
  titlePanel("Adapt GMCP"),
  sidebarLayout(
    sidebarPanel(
      includeCSS("www/styles.css"),
      selectInput("method", "Test Method", choices = c("CombPValue" = "CombPValue", "CER" = "CER"), selected = "CER"),
      numericInput("sampleSize", "Total Sample Size", value = 324),
      numericInput("alpha", "Design Alpha", value = 0.025),
      numericInput("nArms", "Number of Arms", value = 3),
      numericInput("nEps", "Number of Endpoints", value = 2),
      selectInput("testStatCont", "Test Statistics(Cont. Endpoints)", choices = c("z" = "z", "t-equal" = "t-equal", "t-unequal" = "t-unequal"), selected = 't-unequal'),
      selectInput("testStatBin", "Test Statistics(Binary Endpoints)", choices = c("Pooled" = "Pooled", "UnPooled" = "UnPooled"), selected = 'UnPooled'),
      selectInput("fwerControl", "FWER Control Method", choices = c("CombinationTest" = "CombinationTest", "None" = "None"), selected = 'None'),
      h6(tags$b("Endpoint Type(Continuous/Binary")),
      rhandsontableUI("epTypes"),
      h5(tags$b("Arm-wise mean for each Continuous Endpoint")),
      rhandsontableUI("armTable"),
      h5(tags$b("Arm-wise standard deviation for each Continuous Endpoint")),
      rhandsontableUI("sdTable"),
      h5(tags$b("Arm-wise Proportion for each Binary Endpoint")),
      rhandsontableUI("propTable"),
      h5(tags$b("Arm-wise planned allocation ratio")),
      rhandsontableUI("ARTable"),
      h5(tags$b("Correlation between endpoints")),
      matrixInputUI("matrix1"),
      div(class = "highlighted-note",
          "Note: The ordering logic for hypotheses arranges them by iterating
          through each endpoint for every treatment before moving to the next endpoint.
          For example, Hypotheses will follow the following order for a case of 2 Treatment with 3 Endpoints -
          H1 = (Treatment1 vs Control for EP1), H2 = (Treatment2 vs Control for EP1),
          H3 = (Treatment1 vs Control for EP2), H4 = (Treatment2 vs Control for EP2),
          H5 = (Treatment1 vs Control for EP3), H6 = (Treatment2 vs Control for EP3)"),
      h5(tags$b("Initial hypothesis weights")),
      rhandsontableUI("HypWeights"),
      h5(tags$b("Transition matrix")),
      transitionMatrixInputUI("transitionMatrix"),
      uiOutput("testTypeSelect"),
      numericInput("nLooks", "Number of Looks", value = 1),
      h5(tags$b("Information Fraction")),
      IADataInputUI("IAInput"),
      selectInput("typeOfDesign", "Boundary Type",
                  choices = list("O'Brien & Fleming" = 'OF',
                                 "Pocock" = 'P',
                                 "Wang & Tsiatis Delta class" = 'WT',
                                 "Pampallona & Tsiatis" = 'PT',
                                 "Haybittle & Peto" = 'HP',
                                 "Optimum design within Wang & Tsiatis class" = 'WToptimum',
                                 "O'Brien & Fleming type alpha spending" = 'asOF',
                                 "Pocock type alpha spending" = 'asP',
                                 "Kim & DeMets alpha spending" = 'asKD',
                                 "Hwang, Shi & DeCani alpha spending" = 'asHSD',
                                 "No early efficacy stop" = 'noEarlyEfficacy'),
                  selected = "asOF"),
      checkboxInput("CommonStdDev", "Common Standard Deviation(efficacy boundary computation)", value = FALSE),
      checkboxInput("MultipleWinners", "Multiple Winners", value = TRUE),
      checkboxInput("selection", "Enable Treatment Selection", FALSE),
      uiOutput("SelectionInputs"),
      selectInput("implicitSSR", "Implicit SSR",
                  choices = c("Re-allocate samples only from de-selected arms to available arms" = "Selection",
                              "Allocate all the planned samples(for the look) to the available arms" = "All",
                              "No Re-allocation" = "None"),
                  selected = "Selection"),
      numericInput("nSimulation", "Number of Simulations", value = 1000),
      numericInput("seed", "Simulation Seed", value = 100),
      checkboxInput("summaryStat", "Print Summary Statistics File", value = FALSE),
      checkboxInput("plotGraphs", "Plot Initial Graphs", value = TRUE),
      checkboxInput("parallel", "Parallel Computations", value = TRUE),
      checkboxInput("saveFile", "Do you want to save inputs to a CSV?", value = FALSE),
      conditionalPanel(
        condition = "input.saveFile == true",
        textInput("csvDir", "Enter Directory to Save CSV", ""),
        textInput("fileName", "Enter file name", ""),
      ),
      actionButton("simulate", "Simulate", class = "btn-primary"),
    ),
    mainPanel(
      # verbatimTextOutput("listOutput"),
      verbatimTextOutput("result")
    )
  )
)

server <- function(input, output, session) {

  initialData_EPs <- reactive({
    numRows <- input$nEps
    numCols <- 1
    colNames <- "Endpoint Type"

    df <- matrix(rep('Continuous',numRows), nrow = numRows, ncol = numCols)
    df <- as.data.frame(df)
    names(df) <- colNames

    # Pre-fill the first column ("Control") with EP values
    rownames(df) <- paste0("EP", seq_len(nrow(df)))

    df
  })

  initialData_EPsxArms <- reactive({
    numRows <- input$nEps
    numCols <- input$nArms
    colNames <- c("Control", paste0("Treatment", seq_len(numCols - 1)))

    # Create a matrix with 0s for numeric data and fill the first column after converting to a data frame
    df <- matrix(0, nrow = numRows, ncol = numCols)
    df <- as.data.frame(df)
    names(df) <- colNames

    # Pre-fill the first column ("Control") with EP values
    rownames(df) <- paste0("EP", seq_len(nrow(df)))

    # Ensure the rest of the columns are numeric
    for (i in 2:(numCols)) {
      df[[i]] <- as.numeric(df[[i]])
    }

    df
  })

  initialData_Arms <- reactive({
    numRows <- input$nArms
    numCols <- 1
    colNames <- c("AllocationRatio")
    rowNames <- c("Control", paste0("Treatment", seq_len(numRows - 1)))
    df <- matrix(1, nrow = numRows, ncol = numCols )
    colnames(df) <- colNames
    rownames(df) <- rowNames
    for (i in seq_len(1)) {
      df[[i]] <- as.numeric(df[[i]])
    }
    df
  })

  initialData_Hyp <- reactive({
    numRows <- input$nEps * (input$nArms - 1)
    numCols <- 1
    colNames <- c("Weights")
    rowNames <- c(paste0("H", seq_len(numRows)))
    df <- matrix(1, nrow = numRows, ncol = numCols )
    colnames(df) <- colNames
    rownames(df) <- rowNames
    for (i in seq_len(1)) {
      df[[i]] <- as.numeric(df[[i]])
    }
    df
  })
  rhandsontableServer("epTypes", initialData = initialData_EPs)
  rhandsontableServer("armTable", initialData = initialData_EPsxArms)
  rhandsontableServer("sdTable", initialData = initialData_EPsxArms)
  rhandsontableServer("propTable", initialData = initialData_EPsxArms)
  rhandsontableServer("ARTable", initialData = initialData_Arms)
  callModule(matrixInput, "matrix1", dimension = reactive({input$nEps}), simulateTrigger = reactive({input$simulate}))

  rhandsontableServer("HypWeights", initialData = initialData_Hyp)
  callModule(transitionMatrixInput, "transitionMatrix", dimension = reactive({input$nEps * (input$nArms - 1)}))

  output$testTypeSelect <- renderUI({
    if (input$method == "CombPValue") {
      selectInput("testType",
                  "Testing Procedure",
                  choices = c('Bonf', 'Sidak', 'Simes', 'Dunnett', 'Partly-Parametric'),
                  selected = "Partly-Parametric")
    } else if (input$method == "CER") {
      selectInput("testType",
                  "Testing Procedure",
                  choices = c('Parametric', 'Non-Parametric', 'Partly-Parametric'),
                  selected = "Partly-Parametric")
    }
  })
  IADataInputServer("IAInput", reactive(input$nLooks))

  # dropdown creator for Select Endpoint
  endpointChoices <- reactive({
    num <- as.integer(input$nEps)
    choices <- as.character(1:num)
    names(choices) <- choices
    choices <- c(choices, "Overall" = "overall")
    choices
  })
  # Selection related inputs condition on select checkbox input
  output$SelectionInputs <- renderUI({
    if(input$selection) {
      list(
        numericInput("selectionLook", "Selection based on Look", value = 1),

        selectInput("selectionEndPoint", "Selection Endpoint",
                    choices = endpointChoices(),
                    selected = "1"),

        selectInput("selectionScale", "Selection Scale",
                    choices = list("Delta" = "delta",
                                   "Test Statistics" = "teststat",
                                   "Standard Error of the test stat" = "stderror",
                                   "p-value (un-adj)" = "pvalue"),
                    selected = "pvalue"),

        selectInput("selectionCriterion", "Selection Criteria",
                    choices = list("Best r" = "best",
                                   "Threshold for selection" = "threshold",
                                   "Epsilon neighborhood" = "epsilon"),
                    selected = "best"),

        numericInput("selectionParameter", "Selection Parameter for the Selection Criteria", value = 1),

        checkboxInput("keepAssociatedEps", "Keep associated Hypothesis", TRUE)
      )
    }
  })
  observeEvent(input$simulate, {
    # Assuming 'table1-armsTable' is the input ID generated by the module for the rhandsontable
    EpType <- processData_tabularInput(input[['epTypes-tabularInput']])
    Arms.Mean <- processData_tabularInput(input[['armTable-tabularInput']])
    Arms.std.dev <- processData_tabularInput(input[['sdTable-tabularInput']])
    Arms.Prop <- processData_tabularInput(input[['propTable-tabularInput']])
    Arms.alloc.ratio <- unlist(processData_tabularInput(input[['ARTable-tabularInput']]))
    EP.Corr <- processData_corrMatrix(input[['matrix1-corrMatrix']])
    WI <- unlist(processData_tabularInput(input[['HypWeights-tabularInput']]))
    G <- processData_transitionMatrix(input[['transitionMatrix-transitionMatrix']])
    info_frac <- unlist(processData_IAInput(input[['IAInput-dataInput']]))

      alpha = as.numeric(input$alpha)
      SampleSize = as.numeric(input$sampleSize)
      nArms = as.numeric(input$nArms)
      nEps = as.numeric(input$nEps)
      lEpType = EpType
      TestStatCont = input$testStatCont
      TestStatBin = input$testStatBin
      FWERControl = input$fwerControl
      Arms.Mean = Arms.Mean
      Arms.std.dev = Arms.std.dev
      Arms.Prop = Arms.Prop
      Arms.alloc.ratio = Arms.alloc.ratio
      EP.Corr = EP.Corr
      WI = WI
      G = G
      test.type = input$testType
      info_frac = info_frac
      typeOfDesign = input$typeOfDesign
      MultipleWinners = input$MultipleWinners
      CommonStdDev = input$CommonStdDev
      Selection = input$selection
      SelectionLook = as.numeric(input$selectionLook)
      SelectEndPoint = as.numeric(input$selectionEndPoint)
      SelectionScale = input$selectionScale
      SelectionCriterion = input$selectionCriterion
      SelectionParmeter = as.numeric(input$selectionParameter)
      KeepAssosiatedEps = input$keepAssociatedEps
      ImplicitSSR = input$implicitSSR
      nSimulation = as.numeric(input$nSimulation)
      Seed = as.numeric(input$seed)
      SummaryStat = input$summaryStat
      Method = input$method
      plotGraphs = input$plotGraphs
      Parallel = input$parallel

      output$listOutput <- renderPrint({
        str(alpha)
        str(SampleSize)
        str(nArms)
        str(nEps)
        str(TestStat)
        str(FWERControl)
        str(Arms.Mean)
        str(Arms.std.dev)
        str(Arms.alloc.ratio)
        str(EP.Corr)
        str(WI)
        str(G)
        str(test.type)
        str(info_frac)
        str(typeOfDesign)
        str(MultipleWinners)
        str(Selection)
        str(SelectionLook)
        str(SelectEndPoint)
        str(SelectionScale)
        str(SelectionCriterion)
        str(SelectionParmeter)
        str(KeepAssosiatedEps)
        str(ImplicitSSR)
        str(nSimulation)
        str(Seed)
        str(SummaryStat)
        str(Method)
        str(plotGraphs)
        str(Parallel)
      })
      if(input$saveFile){
        dfInputData <- data.frame(Method = Method,
                                  SampleSize = SampleSize,
                                  alpha = alpha,
                                  TestStatCont = TestStatCont,
                                  CommonStdDev = CommonStdDev,
                                  TestStatBin = TestStatBin,
                                  FWERControl = FWERControl,
                                  nArms = nArms,
                                  nEps = nEps,
                                  lEpType = convertToString(lEpType),
                                  Arms.Mean = convertToString(Arms.Mean),
                                  Arms.std.dev = convertToString(Arms.std.dev),
                                  Arms.Prop = convertToString(Arms.Prop),
                                  Arms.alloc.ratio = convertToString(Arms.alloc.ratio),
                                  EP.Corr = convertToString(EP.Corr),
                                  WI = convertToString(WI),
                                  G = convertToString(G),
                                  test.type = test.type,
                                  info_frac = convertToString(info_frac),
                                  typeOfDesign = typeOfDesign,
                                  MultipleWinners = MultipleWinners,
                                  Selection = Selection,
                                  SelectionLook = ifelse(is.null(SelectionLook),"", SelectionLook),
                                  SelectEndPoint = ifelse(is.null(SelectEndPoint),"", SelectEndPoint),
                                  SelectionScale = ifelse(is.null(SelectionScale),"", SelectionScale),
                                  SelectionCriterion = ifelse(is.null(SelectionCriterion),"", SelectionCriterion),
                                  SelectionParmeter = ifelse(is.null(SelectionParmeter),"", SelectionParmeter),
                                  KeepAssosiatedEps = ifelse(is.null(KeepAssosiatedEps),"", KeepAssosiatedEps),
                                  ImplicitSSR = ImplicitSSR,
                                  nSimulation = nSimulation,
                                  Seed = Seed,
                                  SummaryStat = SummaryStat,
                                  plotGraphs = plotGraphs,
                                  Parallel = Parallel)
        dfInputData <- dfInputData %>%
          mutate(ModelID = row_number()) %>%
          select(ModelID, everything())

        csvFileName <- paste0(input$csvDir, "\\", input$fileName,".csv")

        if (file.exists(csvFileName)) {
          existingData <- read.csv(csvFileName)
          combinedData <- rbind(existingData, dfInputData) %>%
            mutate(ModelID = row_number()) %>%
            select(ModelID, everything())
          write.csv(combinedData, csvFileName, row.names = FALSE)
        } else {
          write.csv(dfInputData, csvFileName, row.names = FALSE)
        }
      }

      future({
        result <- simMAMSMEP(alpha = alpha, SampleSize = SampleSize, nArms = nArms, nEps = nEps,lEpType = lEpType,
                   TestStatCont = TestStatCont, TestStatBin = TestStatBin,Arms.Prop = Arms.Prop,CommonStdDev = CommonStdDev,
                   FWERControl = FWERControl, Arms.Mean = Arms.Mean, Arms.std.dev = Arms.std.dev, Arms.alloc.ratio = Arms.alloc.ratio,
                   EP.Corr = EP.Corr,WI = WI, G = G, test.type = test.type, info_frac = info_frac,
                   typeOfDesign = typeOfDesign, MultipleWinners = MultipleWinners,
                   Selection = Selection, SelectionLook = SelectionLook, SelectEndPoint = SelectEndPoint,
                   SelectionScale = SelectionScale, SelectionCriterion = SelectionCriterion,
                   SelectionParmeter = SelectionParmeter, KeepAssosiatedEps = KeepAssosiatedEps,
                   ImplicitSSR = ImplicitSSR, nSimulation = nSimulation, Seed = Seed, SummaryStat = SummaryStat,
                   Method = Method, plotGraphs = plotGraphs, Parallel = Parallel)
        # paste(result)
      }, seed = TRUE) %>%
        # Use then() to handle the promise when it's resolved
        then(function(result) {
          output$result <- renderPrint({
            print(result)
          })
        })
  })
}

shinyApp(ui = ui, server = server)
