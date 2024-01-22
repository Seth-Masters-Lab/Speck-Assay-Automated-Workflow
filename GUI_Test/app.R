library(shiny)
library(flowCore)
rm(list = ls())


files <- list.dirs(path = "Data", recursive = TRUE)

ui <- fluidPage(
  titlePanel("AutoSpeck"),
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = 'data', label = "Select Dataset", choices = files),
      actionButton(inputId = 'confirmBtn', label = "Confirm Selection")
    ),
    mainPanel(
      tableOutput("channels")
    )
  )
)

server <- function(input, output) {
  selectedData <- eventReactive(input$confirmBtn, {
    fcsfiles <- list.files(path = input$data, pattern = "\\.fcs$", ignore.case = TRUE)
    fs <- read.flowSet(fcsfiles[1], input$data, alter.names = T)
    channelValues <- (fs[[1]]@parameters@data[1])
  })
  
  output$channels <- renderTable({
    # Extract elements from the named list
    channelList <- selectedData()
    
    # Create a data frame from the named list
    channelDF <- data.frame(Channels = unlist(channelList), row.names = NULL)
    
    return(channelDF)
  })
  
}


shinyApp(ui = ui, server = server)
