#https://mastering-shiny.org/basic-app.html
library(shiny)
library(Seurat)
library(ggplot2)
load("~/samples/single_cell_samples/single_cell_rnaseq_analysis/samples/SRR8526547/SRR8526547.RData")

reductions <- c("umap","tsne","pca")
ui <- fluidPage(
  sidebarLayout(
  sidebarPanel(
    titlePanel("Seurat two dimensions representation"),
    selectInput("reduction", "Choose a Reduction", choices = reductions),
    sliderInput("x", label = "Resolution:", min = 0.1, max = 1, value = 0.1, step = 0.1)),
  
  mainPanel(
    plotOutput("umap")
)))

server <- function(input, output, session) {
  output$umap <- renderPlot({
    res <- paste0("RNA_snn_res.",input$x)
    Idents(sample) <- res
    DimPlot(sample, reduction=input$reduction, label=TRUE)+labs(title=paste0("resolution: ",res))
  }, width=600, height=600)
}

shinyApp(ui, server)
