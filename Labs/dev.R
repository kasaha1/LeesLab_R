
#====================================functions
## To install Packages-------------
instPak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#------------- Packages ----
packages <- c("ggplot2", "dplyr", "reshape2","data.table","colorspace","dendextend","amap","gplots","corrplot","pca3d","magick","readr")
instPak (packages) 
#-----------------------------


## To install Packages-------------Bioclite

instPak_bioc <- function(pkg_b){
  new.pkg <- pkg_b[!(pkg_b %in% installed.packages()[, "Package"])]
  if (length(new.pkg)){
    source("https://bioconductor.org/biocLite.R")
    biocLite(suppressUpdates=TRUE,suppressAutoUpdate=FALSE,ask=FALSE)
    biocLite(pkg=pkg_b,suppressUpdates=TRUE,suppressAutoUpdate=FALSE,ask=FALSE)
  }
  sapply(pkg_b, require, character.only = TRUE)
}

#------------- Bioc_Packages ----
packages_bioc <- c("ctc")
instPak_bioc (packages_bioc)
#-----------------------------

# Transfomr_NA_to_Median
transform_na_to_median <- function(x) {
  raw.data <- x[-1] %>% as.matrix()
  for (i in c(1:nrow(x))){
    temp.row <- raw.data[i,]
    median.temp <- median(temp.row,na.rm = T)
    raw.data[i,is.na(raw.data[i,])] <- median.temp
  }
  res <- cbind(x[c(1)],raw.data)
  return (res)
}
## Duplicated value removal by SD ---------------
duplicateRemoverbySD <- function(x){
  matrix_data <- as.matrix(x[,-c(1)])
  sd <- apply(matrix_data,1,sd)
  order_num <- seq(1:nrow(x))
  transformed <- cbind(order_num,sd,x)
  name_list <- colnames(transformed)
  colnames(transformed) <- paste0("var_",seq(1:ncol(transformed)))
  colnames(transformed)[1:3] <- c("order_num","sd","grouped")
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

# Gene filtering

geneFilterBySD <- function(x, sdValue =2){
  raw.data <- as.matrix(x[-1])
  sd.filter <- apply(raw.data,1,sd)
  sd_merged <- cbind(sd.filter,x) %>% filter(sd.filter>sdValue)
  print(paste(nrow(sd_merged),"passed out of",nrow(x),"-- FilterBySD"))
  return(sd_merged[-1])
}


geneFilterByABS <- function(x, AtLeastObservation=1, absVal=2){
  raw.data <- as.matrix(x[-1])
  y <- abs(raw.data) %>% apply(c(1,2),function(k){ifelse(k>absVal,TRUE,FALSE)}) %>%  apply(1,sum)
  y_1 <- ifelse(y>=AtLeastObservation,TRUE,FALSE)
  print(paste(nrow(x[y_1,]),"passed out of",nrow(x),"-- FilterByABS"))
  return(x[y_1,])
}
# reorder dend-order
reorder.dend <- function(x) {
  i <- 1
  res.table<- matrix(c(1),ncol = 1,nrow = 1)
  for (j in c(2:(nrow(x)))) {
    if (x[j,2]!=x[(j-1),2]) {i <- i+1}
    temp <- matrix(c(i),nrow = 1,ncol = 1)
    res.table <- rbind(res.table,temp)
  }
  res.table <- as.data.frame(res.table)
  return(cbind(x,res.table))
  # return(cbind(x,res.table))
}

# variables define
library(shiny)
options(shiny.maxRequestSize = 3000 * 1024 ^ 2)

distanceOption <- "pearson"
linkOption <- "average"
tree_cut_r <- 3

status <<- "Please, upload File. <br/> It will be median centered and its NA values will be replaced by median value"
# # # # #


server <- function(input, output, session) {
  # reactive functions 
  datainFile <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    
    f <- fread(inFile$datapath,sep = input$sep, na.strings="NA") %>% as.data.frame() %>% transform_na_to_median %>% duplicateRemoverbySD %>% geneMedianCentering
    print("uploded")
    status <<- "File upload complete. The next step, Data filter, is ready."
    return(f)
  })
  
  reactivefiltering <- reactive({
    data.raw <- datainFile()
    if (is.null(data.raw))
      return(NULL)
    res.filter <- data.raw %>% geneFilterByABS(absVal = input$ABS_value, AtLeastObservation = input$atLeastObservation_val) %>% geneFilterBySD(sdValue = input$SD_val)
    status1 <<- paste("Filtering done. ",nrow(res.filter),"passed out of",nrow(data.raw),"-- FilterByABS")
    status <<- NULL
    return(res.filter)
  })
  
  reactiveHorizontalCluster.hc_dataClu <- reactive({
    d_dataClu <- Dist(raw.data.matrix, method = distanceOption)
    hc_dataClu <- hclust(d_dataClu, method = linkOption)
  })
  reactiveHorizontalCluster.hc_dataClu_r <- reactive({
    d_dataClu_r <- Dist(raw.data.matrix_r, method = distanceOption)
    hc_dataClu_r <- hclust(d_dataClu_r, method = linkOption)
  })
  
  
  # event reactive
  doDataFiltering <- eventReactive(input$doFilter,{
    res.filter <- reactivefiltering()
    
    return(res.filter)
  })
  
  doClustering <- eventReactive(input$doClustering,{
   
    # clustering start
    
    cluster.data <- matrixTranspositionXY(reactivefiltering())
    ID.data <- t(cluster.data[1])
    raw.data.matrix <<- as.matrix(cluster.data[-1])
    
    hc_dataClu <- reactiveHorizontalCluster.hc_dataClu()
    dend <- as.dendrogram(hc_dataClu)
    dend <- color_branches(dend, k= input$tree_cut)
    labels(dend) <- ID.data[order.dendrogram(dend)]
    
    # clutering group to TXT file
    clusterCut <- cutree(hc_dataClu,input$tree_cut)
    ordered_ID <- ID.data[order.dendrogram(dend)] %>% as.data.frame()
    cutting_1 <-  cbind(cluster.data[1],clusterCut)
    file.cutting.sample <<- inner_join(ordered_ID,cutting_1,by=c("."="sample"))
    file.cutting.sample <<- reorder.dend(file.cutting.sample)[-2]
    colnames(file.cutting.sample) <<- c("sample","group")

    ## vertical clustering
    cluster.data_r <- reactivefiltering()
    ID.data_r <- t(cluster.data_r[1])
    raw.data.matrix_r <<- as.matrix(cluster.data_r[-1])
    
    hc_dataClu_r <- reactiveHorizontalCluster.hc_dataClu_r()
   
    dend_r <- as.dendrogram(hc_dataClu_r)

    # print gene
    clusterCut_r <- cutree(hc_dataClu_r,tree_cut_r)
    ordered_ID_r <- ID.data_r[order.dendrogram(dend_r)] %>% as.data.frame()
    cutting_1_r <-  cbind(cluster.data_r[1],clusterCut_r)
    temp_name <- colnames(cutting_1_r)[1]
    file.cutting.gene <<- inner_join(ordered_ID_r,cutting_1_r,by=c("."=temp_name))
    
    
    # heatmap
    
    data.clustering_heatmap <- reactivefiltering()[-1]
    row.names(data.clustering_heatmap) <- t(reactivefiltering()[1])
    # strip_colors <- rainbow_hcl(length(levels(groupData$acronym)))[as.numeric(groupData$acronym)]
    
    
    data_melted <- melt(data.clustering_heatmap,id.vars = 0)
    v_1 <- quantile(data_melted$value,probs = 0.10)
    v_2 <- quantile(data_melted$value,probs = 0.90)
    colors_break = unique(c(seq(min(data_melted$value),v_1,length=100),seq(v_1,v_2,length=100),seq(v_2,max(data_melted$value),length=100)))
    my_palette <- colorRampPalette(c("green", "black", "red"))(n = 297)
    
    heatmap.2(as.matrix(data.clustering_heatmap),
              # main = "Heatmap for the data set",
              # srtCol = 20,
              dendrogram = "both",
              Rowv = dend_r,
              Colv = dend, # this to make sure the columns are not ordered
              density.info = "density",
              breaks=colors_break,
              col = my_palette,
              trace="none",
              cexRow = 0.3,
              symm=F,symkey=F,symbreaks=T, scale="none"
    )
    
    pdf(file="temp_plot.pdf",width = 15, height = 10,pointsize = 12)
      heatmap.2(as.matrix(data.clustering_heatmap),
      # main = "Heatmap for the data set",
      # srtCol = 20,
      dendrogram = "both",
      Rowv = dend_r,
      Colv = dend, # this to make sure the columns are not ordered
      density.info = "density",
      breaks=colors_break,
      col = my_palette,
      trace="none",
      cexRow = 0.3,
      symm=F,symkey=F,symbreaks=T, scale="none"
    )
    dev.off()
  })
  
  # Render
  output$resultUpload <- renderUI({
    # HTML(doUploadFile())
    datainFile()
    HTML(status)
  })
  output$resultfilter <- renderUI({
    doDataFiltering()
    HTML(status1)
  })
  output$cluterPlot <- renderPlot({
    doClustering()
  })

  
  # downLoad
  output$downloadDataFiltered <- downloadHandler(
    filename = function() {
      "FilteredData.txt"
    },
    content = function(file) {
      
      contents.table <- reactivefiltering()
      write_delim(contents.table, file, delim = "\t",na = "")
    },
    contentType = "text/plain"
  )
  output$downloadclusterGroup <- downloadHandler(
    filename = function() {
      "ClusteredSample.txt"
    },
    content = function(file) {
      
      contents.table <- file.cutting.sample
      write_delim(contents.table, file, delim = "\t",na = "")
    },
    contentType = "text/plain"
  )
  output$downloadDataHeatmap <- downloadHandler(
    filename = function() {
      "heatmap.pdf"
    },
    content <- function(file) {
      file.copy("temp_plot.pdf", file)
    }
  )
  
  
}

ui <- fluidPage(titlePanel("UnSupervised Clustering"),
                sidebarLayout(
                  sidebarPanel(
                    fileInput(
                      'file1',
                      h4('Choose txt File'),
                      accept = c('text/csv',
                                 'text/comma-separated-values,text/plain',
                                 '.csv')
                    ),
                    radioButtons('sep', 'Separator',
                                 c(
                                   Comma = ',',
                                   Semicolon = ';',
                                   Tab = '\t'
                                 ),
                                 '\t'),
                    br(),
                    h4("Filter Data"),
                    numericInput("SD_val", label = "SD (gene vector)//", value = 0),
                    numericInput("atLeastObservation_val", label = "At least", value = 1),
                    numericInput("ABS_value", label = "observations with abs(Val)>=", value = 0),
                    actionButton("doFilter", "Filter Data",class = "btn-primary"),
                    downloadButton('downloadDataFiltered', 'Download Filtered Data'),
                    br(),
                    br(),
                    numericInput("tree_cut", label = "Number of clustering group", value = 3),
                    actionButton("doClustering", "Do Cluster", class = "btn-primary"),
                    downloadButton('downloadclusterGroup', 'Download Clustered'),
                    br(),
                    br(),
                    downloadButton('downloadDataHeatmap', 'Download Heatmap'),
                    br(),
                    br()
                  ),
                  mainPanel(
                    # tableOutput("table_display"),
                    htmlOutput("resultUpload"),
                    htmlOutput("resultfilter"),
                    br(),br(),
                    plotOutput("cluterPlot")
                  )
                ))

runApp(shinyApp(ui = ui, server = server),launch.browser = T)