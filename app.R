# Title: Shiny App for Olfactory Expression Data
# Author: Diya Das
# Date: 2018-08-31
# Last revised: Fri Aug 31 12:21:31 2018

# This script needs XCode Developer Tools on MacOS: 
### To install them, open Terminal and run "xcode-select --install"

packagesNeeded = c("shiny", "NMF", "later")
options(repos = structure(c(CRAN="http://cran.cnr.berkeley.edu/")))

new.pkg <- packagesNeeded[!(packagesNeeded %in% installed.packages()[, "Package"])]
if (length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
#based on https://gist.github.com/stevenworthington/3178163

if(!file.exists("oeHBCregenWT_countsMatrix.txt")) unzip("oeHBCregenWT_countsMatrix.txt.zip")

library(shiny)
library(NMF)
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
) #clusterExperiment's bigPalette, avoiding loading due to time
colorscale <- c("#000000", "#00003B", "#000077", "#00009E", "#0000C2", "#0826CD", 
                "#135CCD", "#247AB5", "#39848B", "#378B5C", "#168B26", "#149306", 
                "#5CB21E", "#A0D02E", "#CFE717", "#FFFF00") #clusterExperiment's seqPal5

col.pal <- c("#1B9E77", "cyan", "#E7298A", "darkolivegreen3", "darkorange3",
             "#CCCCCC", "#6A3D9A", "#FCCDE5", "cornflowerblue", "#FFED6F",
             "#FF7F00")
cole <- c("#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594", "#FEE5D9", "#FCBBA1", "#FC9272", "#FB6A4A", "#DE2D26", "#A50F15")

cts <- read.table("oeHBCregenWT_countsMatrix.txt")
clus.labelsdf <- read.table("oeHBCregenWT_clusterLabels.txt")
batchexpt <- read.table("oeHBCregenWT_batchexpt.txt")
# Define UI for application
ui <- fluidPage(
  
  # Application title
  titlePanel("Olfactory gene expression"),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      h3("Select genes here"),
      # Input: Select genes based on data ----
      htmlOutput("multigene_selector"),
      tags$hr(),
      
      # Input: Select a file ----
      fileInput("reflist", "- OR - Choose Reference Gene List",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      checkboxInput("refheader", "Header", FALSE),
      p("A reference gene list will not plot if genes are selected in the box above."),
      tags$hr(),

      htmlOutput("gene_selector"),
      
      h3("Upload custom data"),
      p("Both a counts matrix and cluster labels are required to plot custom datasets."),
      # Input: Select expression file ----
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
      tags$hr(),
      
      # Input: Select a cluster labels file ----
      p("Cluster labels should be uploaded as a two-column text file with sample names in column 1."),
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
      tags$hr(),
      
      # Input: Select a batch and experimental condition file ----
      p("Batch and/or experimental condition should be uploaded in a text file the following headers: sample, batch, expt"),
      fileInput("batchexptdata", "Choose Batch and/or Experiment File",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      radioButtons("batchexptdatasep", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t"),
                   selected = "\t")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      htmlOutput("selected_data"),
      plotOutput("heatmap"),
      plotOutput("lineplot")
    )
  )
)

# Define server logic to read selected file ----
server <- function(input, output) {
  output$selected_data <- renderText({ 
    if (is.null(input$expdata) | is.null(input$clusdata)){
      paste("<h4>Plotting Olfactory Epithelium Regeneration Data</h4>")
    } else {
      paste("<h4>Plotting Uploaded Data</h4>")
    }
  }) 

  output$multigene_selector = renderUI({
    data_available = rownames(cts)
    if (!is.null(input$expdata)){
      cts <- read.table(input$expdata$datapath,
                        header = input$expdataheader,
                        sep = input$expdatasep)
      data_available = rownames(cts)
    }
    selectInput(inputId = "multigene", #name of input
                label = "Genes to Plot in Heatmap", #label displayed in ui
                choices = sort(unique(data_available)),
                multiple = TRUE)
    
  })
  output$gene_selector = renderUI({
    data_available = rownames(cts)
    if (!is.null(input$expdata)){
      cts <- read.table(input$expdata$datapath,
                        header = input$expdataheader,
                        sep = input$expdatasep)
      data_available = rownames(cts)
    }
    selectInput(inputId = "gene", #name of input
                label = "Single Gene to Plot // Doesn't do anything at the moment", #label displayed in ui
                choices = sort(unique(data_available)))
    
  })
  output$heatmap <- renderPlot({
    get_reflist <- reactive({
      if (!is.null(input$multigene)){
        return(input$multigene)
      } else if (!is.null(input$reflist$datapath)) {
        return(read.table(input$reflist$datapath, header = input$refheader))
      } else{
        return(NA)
      }
    })
    
    validate(need(!is.na(get_reflist()), "Can't plot a heatmap without genes! Select a gene from the dropdown OR upload a reference gene list.\n"))
    if (is.null(input$expdata) | is.null(input$clusdata)){

      colb <- bigPalette
      clus.labels <- clus.labelsdf[,2]
      names(clus.labels) <- clus.labelsdf[,1]
      cluster_ord <- c(9,6,2,5,7,11,12,3,8,4,1)
      clus.labels <- clus.labels[order(match(clus.labels,cluster_ord))]
      annot <- data.frame(cluster = factor(clus.labels), batchexpt)
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
      
      if (!is.null(input$batchexptdata)){
        batchexpt <- read.table(input$batchexptdata$datapath,
                                header = TRUE,
                                sep = input$batchexptdatasep,
                                row.names = 1)
        batchexpt <- batchexpt[names(clus.labels),]
        annot <- data.frame(cluster = factor(clus.labels), batchexpt)
      } else {
        annot <- data.frame(cluster = factor(clus.labels))
      }
    }
    breakv <- c(min(cts),
                seq(0, quantile(cts[cts > 0], .99, na.rm = TRUE), length = 50),
                max(cts))
    ph <- aheatmap(cts[intersect(as.character(unlist(get_reflist())), rownames(cts)), names(clus.labels)], Rowv = NA, Colv = NA, breaks = breakv, annCol = annot, annColors = list(cluster = col.pal, expt = cole, batch = colb), color = colorscale, annLegend = TRUE)
    return(ph)
  })
  
  output$lineplot <- renderPlot({
    req(input$gene)
    if (is.null(input$expdata) | is.null(input$clusdata)){
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
      clus.labels <- clus.labelsdf[,2]
      names(clus.labels) <- clus.labelsdf[,1]
      clus.labels <- sort(clus.labels)
      
      if (!is.null(input$batchexptdata)){
        batchexpt <- read.table(input$batchexptdata$datapath,
                                header = TRUE,
                                sep = input$batchexptdatasep,
                                row.names = 1)
        batchexpt <- batchexpt[names(clus.labels),]
      } else{
        batchexpt <- NULL
      }
      
    }
    par(mfrow=c(3,1), mar=c(4,4,2,1))
    plot(t(cts[input$gene, names(clus.labels)]), col = col.pal[clus.labels], pch = 19, xlab = "By Cluster", ylab = input$gene)
    plot(t(cts[input$gene, names(clus.labels)]), col = cole[batchexpt$expt], pch = 19, xlab = "By Experimental Condition", ylab = input$gene)
    plot(t(cts[input$gene, names(clus.labels)]), col = bigPalette[batchexpt$batch], pch = 19, xlab = "By Batch", ylab = input$gene)
  })
}


# Run the app ----
shinyApp(ui, server)
