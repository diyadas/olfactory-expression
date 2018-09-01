# Title: Shiny App for Olfactory Expression Data
# Author: Diya Das
# Date: 2018-08-31
# Last revised: Fri Aug 31 12:21:31 2018

# This script needs XCode Developer Tools on MacOS: 
### To install them, open Terminal and run "xcode-select --install"

packagesNeeded = c("shiny", "NMF", "later")
options(repos = structure(c(CRAN="http://cran.cnr.berkeley.edu/")))

new.pkg <- packagesNeeded[!(packagesNeeded %in% 
                              installed.packages()[, "Package"])]
if (length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
#based on https://gist.github.com/stevenworthington/3178163

if(!file.exists("oeHBCregenWT_countsMatrix.txt")) unzip("oeHBCregenWT_countsMatrix.txt.zip")

library(shiny)
library(NMF)
options(shiny.maxRequestSize=1024^3) 
colb <- c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A", 
                "#B15928", "#A6CEE3", "#bd18ea", "cyan", "#B2DF8A", "#FB9A99", 
                "deeppink4", "#00B3FFFF", "#CAB2D6", "#FFFF99", "#05188a", 
                "#CCFF00FF", "cornflowerblue", "#f4cc03", "black", 
                "blueviolet", "#4d0776", "maroon3", "blue", "#E5D8BD", 
                "cadetblue4", "#e5a25a", "lightblue1", "#F781BF", "#FC8D62", 
                "#8DA0CB", "#E78AC3", "green3", "#E7298A", "burlywood3", 
                "#A6D854", "firebrick", "#FFFFCC", "mediumpurple", "#1B9E77", 
                "#FFD92F", "deepskyblue4", "yellow3", "#00FFB2FF", "#FDBF6F", 
                "#FDCDAC", "gold3", "#F4CAE4", "#E6F5C9", "#FF00E6FF", 
                "#7570B3", "goldenrod", "#85848f", "lightpink3", "olivedrab",
                "cadetblue3") #clusterExperiment's bigPalette
colscale <- c("#000000", "#00003B", "#000077", "#00009E", "#0000C2", "#0826CD",
              "#135CCD", "#247AB5", "#39848B", "#378B5C", "#168B26", "#149306",
              "#5CB21E", "#A0D02E", "#CFE717", "#FFFF00") 
                #clusterExperiment's seqPal5
col.pal <- c("cornflowerblue", "#CCCCCC", "cyan", "darkorange3",  "#6A3D9A",
             "#FFED6F", "#FF7F00", "#E7298A", "#FCCDE5", "darkolivegreen3",
             "#1B9E77")
cole <- c("#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594", 
          "#FEE5D9", "#FCBBA1", "#FC9272", "#FB6A4A", "#DE2D26", "#A50F15")

cts <- read.table("oeHBCregenWT_countsMatrix.txt")
clus.labelsdf <- read.table("oeHBCregenWT_clusterLabels.txt", sep = "\t")
clus.labels <- clus.labelsdf[,2]
names(clus.labels) <- clus.labelsdf[,1]
cluster_ord <- c("activated HBC 1", "activated HBC 2", "GBC", "INP1/2", "INP3",
                 "iOSN", "mOSN", "SUS", "renewing HBC 1",  "renewing HBC 2",
                 "resting HBC")
clus.labels <- factor(clus.labels, levels = cluster_ord)
clus.labels <- clus.labels[order(match(clus.labels,cluster_ord))]
batchexpt <- read.table("oeHBCregenWT_batchexpt.txt", sep = "\t", 
                        header = TRUE, row.names = 1)
batchexpt$expt <- factor(batchexpt$expt, 
                         levels = c("Uninjured", "24 HPI", "48 HPI",
                                    "96 HPI", "7 DPI", "14 DPI"))
batchexpt <- batchexpt[names(clus.labels),]

# Define UI for application
ui <- fluidPage(
  tags$style(".shiny-input-container { margin-bottom: 0px; margin-top: 0px } 
             #reflist_progress { margin-bottom: 0px } #clusdata_progress { 
             margin-bottom: 0px } #expdata_progress { margin-bottom: 0px } 
             #batchexptdata_progress { margin-bottom: 0px } .checkbox { 
             margin-top: 0px }"),
  
  # Application title
  titlePanel("Olfactory gene expression"),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      h3("Select genes here"),
      # Input: Select genes based on data ----
      htmlOutput("multigene_selector"),
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
      tags$hr(),
      
      h3("Upload custom data"),
      p("Both a counts matrix and cluster labels are required to plot custom 
        datasets."),
      # Input: Select expression file ----
      fileInput("expdata", "Choose Expression Data File",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      radioButtons("expdatasep", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t"),
                   selected = "\t"),
      
      # Input: Select a cluster labels file ----
      p("Cluster labels should be uploaded as a two-column text file 
        with sample names in column 1."),
      fileInput("clusdata", "Choose Cluster Labels File",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      checkboxInput("clusdataheader", "Header", FALSE),
      radioButtons("clusdatasep", "Separator",
                   choices = c(Comma = ",", Semicolon = ";", Tab = "\t"),
                   selected = "\t"),

      # Input: Select a batch and experimental condition file ----
      p("Batch and/or experimental condition should be uploaded in a text 
        file with the following column names: sample, batch, expt"),
      fileInput("batchexptdata", "Choose Batch and/or Experiment File 
                (Optional)",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      radioButtons("batchexptdatasep", "Separator",
                   choices = c(Comma = ",", Semicolon = ";", Tab = "\t"),
                   selected = "\t")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      htmlOutput("selected_data"),
      plotOutput("heatmap", height = "600px"),
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
                        header = TRUE,
                        sep = input$expdatasep)
      data_available = rownames(cts)
    }
    selectInput(inputId = "multigene", 
                label = "Genes to Plot in Heatmap", 
                choices = sort(unique(data_available)),
                multiple = TRUE)
    
  })
  output$gene_selector = renderUI({
    data_available = rownames(cts)
    if (!is.null(input$expdata)){
      cts <- read.table(input$expdata$datapath,
                        header = TRUE,
                        sep = input$expdatasep)
      data_available = rownames(cts)
    }
    selectInput(inputId = "gene", 
                label = "Single Gene to Plot", 
                choices = sort(unique(data_available)))
  })
  output$heatmap <- renderPlot({
    get_reflist <- reactive({
      if (!is.null(input$multigene)) {
        return(input$multigene)
      } else if (!is.null(input$reflist$datapath)) {
        return(read.table(input$reflist$datapath, header = input$refheader))
      } else {
        return(NA)
      }
    })
    validate(need(!is.na(get_reflist()), 
                  paste("Can't plot a heatmap without genes! Select at least", 
                        "one gene OR upload a reference gene list.")))
    
    if (is.null(input$expdata) | is.null(input$clusdata)){
      annot <- data.frame(cluster = factor(clus.labels), 
                          batchexpt[,2:1])
    } else {
      cts <- read.table(input$expdata$datapath,
                        header = TRUE,
                        sep = input$expdatasep)
      clus.labelsdf <- read.table(input$clusdata$datapath,
                                  header = input$clusdataheader,
                                  sep = input$clusdatasep)
      col.pal <- cole <- colb
      clus.labels <- clus.labelsdf[,2]
      names(clus.labels) <- clus.labelsdf[,1]
      clus.labels <- sort(clus.labels)
      if (!is.null(input$batchexptdata)){
        batchexpt <- read.table(input$batchexptdata$datapath,
                                header = TRUE,
                                sep = input$batchexptdatasep,
                                row.names = 1)
        annot <- data.frame(cluster = clus.labels, 
                            batchexpt[names(clus.labels), ])
      } else {
        annot <- data.frame(cluster = clus.labels)
      }
    }
    breakv <- c(min(cts),
                seq(0, quantile(cts[cts > 0], .99, na.rm = TRUE), length = 50),
                max(cts))
    ph <- aheatmap(cts[intersect(as.character(unlist(get_reflist())), 
                                 rownames(cts)), names(clus.labels)], 
                   Rowv = NA, Colv = NA, breaks = breakv, annCol = annot,
                   annColors = list(cluster = col.pal, expt = cole, 
                                    batch = colb), 
                   color = colscale, annLegend = TRUE)
    return(ph)
  })
  
  output$lineplot <- renderPlot({
    req(input$gene)
    if (!is.null(input$expdata) & !is.null(input$clusdata)){
      cts <- read.table(input$expdata$datapath,
                        header = TRUE,
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
      } 
    }
    par(mfrow=c(3,1), mar=c(4,4,2,1))
    plot(t(cts[input$gene, names(clus.labels)]), col = col.pal[clus.labels],
         pch = 19, xlab = "By Cluster", ylab = input$gene)
    plot(t(cts[input$gene, names(clus.labels)]), col = cole[batchexpt$expt],
         pch = 19, xlab = "By Experimental Condition", ylab = input$gene)
    plot(t(cts[input$gene, names(clus.labels)]), col = colb[batchexpt$batch], 
         pch = 19, xlab = "By Batch", ylab = input$gene)
  })
}

# Run the app ----
shinyApp(ui, server)
