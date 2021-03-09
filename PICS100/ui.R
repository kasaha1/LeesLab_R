#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    # Application title
    titlePanel("Temporary PICS100 prediction"),
    
    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            fileInput(
                'geneExprfile',
                h4('Choose Gene expression file(txt/csv)'),
                accept = c('text/csv',
                           'text/comma-separated-values,text/plain',
                           '.csv')
            ),
            radioButtons(
                'standardizationType',
                'Standardization',
                c(
                    'Mendian-centering only' = 'medianCenter',
                    'Median-centering and dividing by SD' = 'devidedBySD'
                ),
                'medianCenter'
            ),
            actionButton("doPrediction", "Prediction", class = "btn-primary"),
            br(),
            br(),
            downloadButton('downloadResults', 'Download Result')
        ),
        
        # Show a plot of the generated distribution
        mainPanel(tabsetPanel(
            type = "tabs",
            tabPanel("Analysis", textOutput("status"),
                     tableOutput("tablesTemp")),
            tabPanel(
                "About PICS100",
                HTML(
                    "<p>While many studies revealed genomic subtypes of hepatocellular carcinoma (HCC), they are not translated to the clinic yet due to lack of consensus. We aim to examine consensus of genomic subtypes and uncover their clinical significance. We integrated 15 previously established genomic signatures for HCC to uncover consensus genomic subtypes.<b> We also developed and validated a robust predictor of consensus subtype with 100 genes (PICS100).<b><p> "
                ),
                img(
                    src = "Fig2.png",
                    width = 500,
                    height = 600
                )
            )
        ))
    )
))
