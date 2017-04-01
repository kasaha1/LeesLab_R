

## To install Packages-------------
instPak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#------------- Packages ----
packages <-
  c("ggplot2",
    "dplyr",
    "reshape2",
    "moonBook",
    "readr",
    "colorspace",
    "shiny")
instPak (packages)
#-----------------------------
if (!("limma" %in% installed.packages()[, "Package"])) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("limma")
}
library(limma)
options(shiny.maxRequestSize = 3000 * 1024 ^ 2)

# # # # #

server <- function(input, output) {
  datainFile <- reactive({
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    read_delim(inFile$datapath,
               delim = input$sep,
               quote = input$quote)
  })
  
  do_it <-  eventReactive(input$action1, {
    data.plot <- datainFile()
    M.nor <- data.plot[-1] %>% as.matrix()
    
    if (input$Flooring_after_Log2) {
      flooringFunction <- function(x) {
        if (x < 1) {
          x <- 1
        }
        return(x)
      }
      M.nor <- apply(M.nor, c(1, 2), flooringFunction)
    }
    if (input$doLog2)
    {
      M.nor <- log2(M.nor)
    }
    if (input$doNoramlization) {
      M.nor <- normalizeBetweenArrays(M.nor, method = "quantile")
    } else{
      M.nor <- M.nor
    }
    return(round(M.nor, 4))
  })
  
  
  output$plotContents1 <- renderPlot({
    if (!is.null(datainFile())) {
      data.plot <- datainFile()
      boxplot(data.plot[-1], main = "Before Do")
      geneName <<- data.plot[1]
    }
  })
  output$plotContents2 <- renderPlot({
    data.plot.do <- do_it()
    boxplot(data.plot.do, main = "After Do")
    genomatrix <<- data.plot.do
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("AfterNormalization", '.txt', sep = '')
    },
    content = function(file) {

      contents.table <- cbind(datainFile()[1], do_it())
      write_delim(contents.table, file, delim = "\t")
    },
    contentType = "text/plain"
  )
}

ui <- fluidPage(titlePanel("Uploading Files"),
                sidebarLayout(
                  sidebarPanel(
                    fileInput(
                      'file1',
                      'Choose txt File',
                      accept = c('text/csv',
                                 'text/comma-separated-values,text/plain',
                                 '.csv')
                    ),
                    tags$hr(),
                    radioButtons('sep', 'Separator',
                                 c(
                                   Comma = ',',
                                   Semicolon = ';',
                                   Tab = '\t'
                                 ),
                                 '\t'),
                    radioButtons(
                      'quote',
                      'Quote',
                      c(
                        None = '',
                        'Double Quote' = '"',
                        'Single Quote' = "'"
                      ),
                      ''
                    ),
                    h4("After file upload......."),
                    checkboxInput("doNoramlization", label = "Noramlization", value = TRUE),
                    checkboxInput("Flooring_after_Log2", label = "Flooring under 1", value = F),
                    checkboxInput("doLog2", label = "Log2", value = F),
                    actionButton("action1", "Do", class = "btn-primary"),
                    downloadButton('downloadData', 'Download')
                  ),
                  mainPanel(
                    plotOutput('plotContents1'),
                    br(),
                    plotOutput('plotContents2')
                  )
                ))
shinyApp(ui = ui, server = server)
