options(repos = BiocManager::repositories())
#setOption("repos")

library(shiny)
library(clusterExperiment)
options(shiny.maxRequestSize=1024^3) 
bigPalette <- c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A", "#B15928", 
                "#A6CEE3", "#bd18ea", "cyan", "#B2DF8A", "#FB9A99", "deeppink4", 
                "#00B3FFFF", "#CAB2D6", "#FFFF99", "#05188a", "#CCFF00FF", "cornflowerblue", 
                "#f4cc03", "black", "blueviolet", "#4d0776", "maroon3", "blue", 
                "#E5D8BD", "cadetblue4", "#e5a25a", "lightblue1", "#F781BF", 
                "#FC8D62", "#8DA0CB", "#E78AC3", "green3", "#E7298A", "burlywood3", 
                "#A6D854", "firebrick", "#FFFFCC", "mediumpurple", "#1B9E77", 
                "#FFD92F", "deepskyblue4", "yellow3", "#00FFB2FF", "#FDBF6F", 
                "#FDCDAC", "gold3", "#F4CAE4", "#E6F5C9", "#FF00E6FF", "#7570B3", 
                "goldenrod", "#85848f", "lightpink3", "olivedrab", "cadetblue3"
) #clusterExperiment's bigPalette, in case clusterExperiment doesn't work


# Define UI for application
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
    checkboxInput("expdataheader", "Header", TRUE),
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
    checkboxInput("clusdataheader", "Header", FALSE),
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
    checkboxInput("refheader", "Header", FALSE),
    
    # Horizontal line ----
    tags$hr(),
    
    # Input: Select genes based on data ----
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
    req(input$reflist)

    reflist <- read.table(input$reflist$datapath, header = input$refheader)
    if (is.null(input$expdata) | is.null(input$clusdata)){
      tmp <- tempfile()
      download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE99251&format=file&file=GSE99251%5FoeHBC%5FWTregen%5FcountsMatrix%2Etxt%2Egz",tmp)
      cts <- read.table(gzfile(tmp))
      unlink(tmp)
      clus.labelsdf <- read.table("https://raw.githubusercontent.com/diyadas/HBC-regen/master/ref/oeHBCregenWT_clusterLabels.txt")
      col.pal <- c("#1B9E77", "cyan", "#E7298A", "darkolivegreen3", "darkorange3", "#CCCCCC", "#6A3D9A", "#FCCDE5", "cornflowerblue", "#FFED6F", "#FF7F00")
      cole <- c("#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594", "#FEE5D9", "#FCBBA1", "#FC9272", "#FB6A4A", "#DE2D26", "#A50F15")
      colb <- bigPalette
      clus.labels <- clus.labelsdf[,2]
      names(clus.labels) <- clus.labelsdf[,1]
      cluster_ord <- c(9,6,2,5,7,11,12,3,8,4,1)
      clus.labels <- clus.labels[order(match(clus.labels,cluster_ord))]
    } else {
      cts <- read.table(input$expdata$datapath,
                        header = input$expdataheader,
                        sep = input$expdatasep)
      clus.labelsdf <- read.table(input$clusdata$datapath,
                                  header = input$clusdataheader,
                                  sep = input$clusdatasep)
      col.pal <- cole <- colb <- bigPalette
      clus.labels <- clus.labelsdf[,2]
      names(clus.labels) <- clus.labelsdf[,1]
      clus.labels <- sort(clus.labels)
    }
    
    
    
    breakv <- c(min(cts),
                seq(0, quantile(cts[cts > 0], .99, na.rm = TRUE), length = 50),
                max(cts))


    #plotHeatmap(cts[intersect(reflist, rownames(cts)), names(clus.labels)], clusterSamples = FALSE, clusterFeatures = FALSE, breaks = breakv, colData = data.frame(cluster = clus.labels, expt = expt, batch = batch), clusterLegend = list(cluster = bigPalette, expt = cole), annLegend = TRUE)

    
    ph <- plotHeatmap(cts[intersect(as.character(unlist(reflist)), rownames(cts)), names(clus.labels)], clusterSamples = FALSE, clusterFeatures = FALSE, breaks = breakv, colData = data.frame(cluster = clus.labels), clusterLegend = list(cluster = col.pal), annLegend = TRUE)
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
