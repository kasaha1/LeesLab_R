ui <- fluidPage(titlePanel(""),
                sidebarLayout(
                  sidebarPanel(
                  
                  ),
                  mainPanel(
                    
                  )
                ))


server <- function(input, output, session) {
 
}

runApp(shinyApp(ui = ui, server = server),launch.browser = T)
