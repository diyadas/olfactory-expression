options(repos = BiocManager::repositories())
#setOption("repos")

library(shiny)
library(clusterExperiment)
options(shiny.maxRequestSize=1024^3) 

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
titlePanel("Olfactory gene expression"),
# Sidebar layout with input and output definitions ----
sidebarLayout(
  # Sidebar panel for inputs ----
  sidebarPanel(
    # Input: Select a file ----
    fileInput("expdata", "Choose Expression Data File",
              multiple = TRUE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),
    
    # Input: Checkbox if file has header ----
    checkboxInput("expdataheader", "Header", TRUE),
    
    # Input: Select separator ----
    radioButtons("expdatasep", "Separator",
                 choices = c(Comma = ",",
                             Semicolon = ";",
                             Tab = "\t"),
                 selected = "\t"),
    
    # Horizontal line ----
    tags$hr(),
    
    # Input: Select a file ----
    fileInput("clusdata", "Choose Cluster Labels File",
              multiple = TRUE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),
    
    # Input: Checkbox if file has header ----
    checkboxInput("clusdataheader", "Header", FALSE),
    
    # Input: Select separator ----
    radioButtons("clusdatasep", "Separator",
                 choices = c(Comma = ",",
                             Semicolon = ";",
                             Tab = "\t"),
                 selected = "\t"),
    
    # Horizontal line ----
    tags$hr(),

    # Input: Select a file ----
    fileInput("reflist", "Choose Reference Gene List",
              multiple = TRUE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),
    
    # Input: Checkbox if file has header ----
    checkboxInput("refheader", "Header", FALSE),
    
    # Horizontal line ----
    tags$hr(),
    
    htmlOutput("gene_selector")
  ),
  
  # Main panel for displaying outputs ----
  mainPanel(
    
    # Output: Data file ----
    tableOutput("contents"),
    tableOutput("clus"),
    plotOutput("heatmap")#,
    #plotOutput("lineplot")
  )
  
)
)

# Define server logic to read selected file ----
server <- function(input, output) {
  
#output$heatmap <- renderPlot(plot(1:3,4:6))
  output$contents <- renderTable({
    req(input$reflist)
    reflist <- read.table(input$reflist$datapath, header = input$refheader)
    return(head(reflist))
  })
  
  output$clus <- renderTable({
    req(input$clusdata)
    clus.labelsdf <- read.table(input$clusdata$datapath,
                                header = input$clusdataheader,
                                sep = input$clusdatasep)
    clus.labels <- clus.labelsdf[,2]
    names(clus.labels) <- clus.labelsdf[,1]
    return(head(clus.labelsdf))
  })

  output$heatmap <- renderPlot({
    req(input$expdata)
    req(input$reflist)
    req(input$clusdata)
    reflist <- read.table(input$reflist$datapath, header = input$refheader)
    
    cts <- read.table(input$expdata$datapath,
                         header = input$expdataheader,
                         sep = input$expdatasep)
    clus.labelsdf <- read.table(input$clusdata$datapath,
                      header = input$clusdataheader,
                      sep = input$clusdatasep)
    clus.labels <- clus.labelsdf[,2]
    names(clus.labels) <- clus.labelsdf[,1]
    breakv <- c(min(cts),
                seq(0, quantile(cts[cts > 0], .99, na.rm = TRUE), length = 50),
                max(cts))


    #plotHeatmap(cts[intersect(reflist, rownames(cts)), names(clus.labels)], clusterSamples = FALSE, clusterFeatures = FALSE, breaks = breakv, colData = data.frame(cluster = clus.labels, expt = expt, batch = batch), clusterLegend = list(cluster = bigPalette, expt = cole), annLegend = TRUE)

    clus.labels <- sort(clus.labels)
    ph <- plotHeatmap(cts[intersect(as.character(unlist(reflist)), rownames(cts)), names(clus.labels)], clusterSamples = FALSE, clusterFeatures = FALSE, breaks = breakv, colData = data.frame(cluster = clus.labels), clusterLegend = list(cluster = bigPalette), annLegend = TRUE)
    return(ph)
  })
  
  output$gene_selector = renderUI({
    req(input$expdata)
    cts <- read.table(input$expdata$datapath,
                      header = input$expdataheader,
                      sep = input$expdatasep)
    data_available = rownames(cts)
    selectInput(inputId = "gene", #name of input
                label = "Gene to Plot", #label displayed in ui
                choices = unique(data_available))
  })
  
  # output$lineplot <- renderPlot({
  #   req(input$expdata)
  #   req(input$gene)
  #   req(input$clusdata)
  #   
  #   cts <- read.table(input$expdata$datapath,
  #                     header = input$expdataheader,
  #                     sep = input$expdatasep)
  #   clus.labelsdf <- read.table(input$clusdata$datapath,
  #                               header = input$clusdataheader,
  #                               sep = input$clusdatasep)
  #   clus.labels <- clus.labelsdf[,2]
  #   names(clus.labels) <- clus.labelsdf[,1]
  #   plot(cts[input$gene,names(clus.labels)],col=bigPalette[clus.labels], pch=19)
  #     lines(lowess(cts[input$gene,],f=0.15,delta=2),lwd=2)
  # })
}


# Run the app ----
shinyApp(ui, server)
