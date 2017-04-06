
## To install Packages-------------
instPak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#------------- Packages ----
packages <- c("readr","dplyr","shiny","data.table")
instPak (packages) 
#-----------------------------

## Duplicated value removal by SD ---------------
duplicateRemoverbySD <- function(x,startCol=1){
  matrix_data <- as.matrix(x[,-c(1:startCol)])
  sd <- apply(matrix_data,1,sd,na.rm=T)
  order_num <- seq(1:nrow(x))
  transformed <- cbind(order_num,sd,x)
  name_list <- colnames(transformed)
  colnames(transformed) <- paste0("var_",seq(1:ncol(transformed)))
  colnames(transformed)[1:2] <- c("order_num","sd")
  colnames(transformed)[(startCol+2)] <- "grouped"
  res <- transformed %>% arrange(desc(sd)) %>% group_by(grouped) %>% filter(row_number()==1) %>% ungroup() %>% arrange(order_num)
  colnames(res) <- name_list
  return(res[c(-1,-2)])
}
## Transpostion XY----------------
matrixTranspositionXY <- function(x, firstColumnName="sample"){
  col_names_1 <- t(x[1])
  raw_data <- t(x[-1])
  colnames(raw_data) <- col_names_1
  raw_data <- as.data.frame(raw_data)
  row_name_1 <- row.names(raw_data)
  raw_data <- cbind(row_name_1,raw_data)
  row.names(raw_data) <- NULL
  colnames(raw_data)[1] <- firstColumnName
  raw_data[,1] <- as.character(raw_data[,1])
  return(raw_data)
}
# gene median centering
geneMedianCentering <- function(x){
  raw.data <- x[-1] %>% as.matrix()
  median.table <- apply(raw.data ,c(1),median,na.rm = T) 
  median_centered <- raw.data-median.table
  return(cbind(x[1],median_centered))
}


options(shiny.maxRequestSize = 3000 * 1024 ^ 2)





ui <-
  fluidPage(
    titlePanel("Extract Gene siganture from the phenotype"),
    sidebarLayout(
      sidebarPanel(
        fileInput(
          'fileRawdata',
          'Choose Gene Expression File',
          accept = c('text/csv',
                     'text/comma-separated-values,text/plain',
                     '.csv')
        ),
        radioButtons('sep1', 'Separator',
                     c(
                       Comma = ',',
                       Semicolon = ';',
                       Tab = '\t'
                     ),
                     '\t'),
        fileInput(
          'fileGroup',
          'Choose Group File',
          accept = c('text/csv',
                     'text/comma-separated-values,text/plain',
                     '.csv')
        ),
        radioButtons('sep2', 'Separator',
                     c(
                       Comma = ',',
                       Semicolon = ';',
                       Tab = '\t'
                     ),
                     '\t'),
        actionButton("actionExtractColumns", "Group extraction", class = "btn-primary"),
        
        selectInput("actionGroup", "Select Group", choices = NULL),
        numericInput("pValue", label = "P value Cutoff", value = 0.001,width = "50%"),
        numericInput("differenceRatio", label = "Group Difference ration", value = 2,width = "50%"),
        
        actionButton("actionExtract", "Extract Signature", class = "btn-primary"),
        br(),
        br(),
        downloadButton('downloadData', 'Download Result')
      ),
      mainPanel(
        htmlOutput("resultFileUpdate"),
        htmlOutput("FilteredGenes")
      )
    )
  )


server <- function(input, output, session) {
  # data reactive
  datainFileRawdata <- reactive({
    inFile <- input$fileRawdata
    req(inFile)
    f <- fread(inFile$datapath,sep = input$sep1, na.strings="NA") %>% as.data.frame()
    print("raw data uploded")
    return(f)
  })
  datainFileGroup <- reactive({
    inFile <- input$fileGroup
    req(inFile)
    f <- fread(inFile$datapath,sep = input$sep2, na.strings="NA") %>% as.data.frame()
    print("group data uploded")
    return(f)
  })
  data.groupRearranged <- reactive({
    raw.groupdata <- datainFileGroup()
    colnames(raw.groupdata)[c(1,2)] <- c("sample","group")
    raw.groupdata$group <- as.character(raw.groupdata$group)
    raw.groupdata$group[raw.groupdata$group!=input$actionGroup] <- "zTheOthers"
    raw.groupdata$group <- as.character(raw.groupdata$group)
    raw.groupdata$group <- as.factor(raw.groupdata$group)
    print("group rearranged")
    return(raw.groupdata)
  })
  
  
  # eventReactive
  select.column <- eventReactive(input$actionExtractColumns, {
    datainFileRawdata()
    fgroup <- datainFileGroup()
    vars <- as.factor(fgroup[,c(2)]) %>% levels
    print("vectors are determinded")
    # Update select input immediately after clicking on the action button. 
    updateSelectInput(session, "actionGroup","Select Group", choices = vars)
    return("File update complete")
  })
  doExtractingSignatutre <- eventReactive(input$actionExtract,{
    data.process <- datainFileRawdata()
    data.process_1<- matrixTranspositionXY(data.process)
    groupData <- data.groupRearranged()
## ERROR 
    a_group <- groupdata %>% filter(group==input$actionGroup)
    b_group <- groupdata %>% filter(group=="zTheOthers")
    print("group divided")
    data.process_1_a <- inner_join(a_group,data.process_1,by=c("sample"="sample"))[-2] %>% matrixTranspositionXY()
    data.process_1_b <- inner_join(b_group,data.process_1,by=c("sample"="sample"))[-2] %>% matrixTranspositionXY()
    res.table <- matrix(nrow = 1,ncol = 2)
    data.process_1_a_m <- data.process_1_a[c(2:ncol(data.process_1_a))] %>% as.matrix()
    data.process_1_b_m <- data.process_1_b[c(2:ncol(data.process_1_b))] %>% as.matrix()
    for (i in 1:nrow(data.process_1_a)) {
      k_1 <- data.process_1_a_m[i,] 
      k_2 <- data.process_1_b_m[i,] 
      average_d <- mean(k_1,na.rm = T)-mean(k_2,na.rm = T)
      t_p <- t.test(k_1,k_2)
      t_p <- t_p$p.value
      res.tem <- matrix(c(average_d,t_p),nrow = 1,ncol = 2)
      res.table <- rbind(res.table,res.tem)
    }
    res.table <- res.table[-1,] %>% as.data.frame()
    res.table <<- cbind(data.process_1_a[1],res.table)
    colnames(res.table) <- c("gene","average_d","t_test_p")
    
    
    
    
    HTML("extraction done")
  })
  # ouput
  output$resultFileUpdate <- renderUI({
    HTML(select.column())
  })
  
  output$FilteredGenes <- renderUI({
    HTML(doExtractingSignatutre())
  })
  
 
}

runApp(shinyApp(ui = ui, server = server),launch.browser = T)