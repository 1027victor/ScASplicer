library(MARVEL)

# Load adjunct packages for selected MARVEL features
# General data processing, plotting
library(ggnewscale)
library(ggrepel)
library(parallel)
library(reshape2)
library(stringr)
library(textclean)

# Dimension reduction analysis
library(factoextra)
library(FactoMineR)

# Modality analysis
library(fitdistrplus)

# Differential splicing analysis
library(kSamples)
library(twosamples)

# Gene ontology analysis
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

# Nonsense-mediated decay (NMD) analysis
library(Biostrings)
library(BSgenome)
library(BSgenome.Hsapiens.NCBI.GRCh38)

# Load adjunct packages for this tutorial
library(data.table)
library(ggplot2)
library(gridExtra)
library(highcharter)
library(shinycssloaders)
library(dplyr)
library(scales)
library(stringr)
library(tidyr) 
library(circlize)                                     
library(grid)
library(graphics)
library(ComplexHeatmap)
# library(openxlsx)

options(future.globals.maxSize = 300* 1024^3)

function_file_lst <- list.files("functions", pattern = ".R$", full.names = TRUE)
for (i in function_file_lst) {
  source(i)
}

module_file_lst <- list.files("module", pattern = ".R$", full.names = TRUE)
for (i in module_file_lst) {
  source(i)
}

shinyServer(function(input, output,session) {
  # Overview of splicing events
  options(shiny.maxRequestSize=300*1024^3)
  
  observeEvent(session$clientData$url_hash, {
    currentHash <- sub("#", "", session$clientData$url_hash)
    if(is.null(input$navBar) || !is.null(currentHash) && currentHash != input$navBar){
      freezeReactiveValue(input, "navBar")
      updateNavbarPage(session, "navBar", selected = currentHash)
    }
    switch(currentHash,
           "event" = {
             overview_mainServer("event")
           },
           "ma" = {
             modal_mainServer("modal")
           },
           "analysis" = {
             Differential_mainServer("diff")
           },
           "go" = {
             go_mainServer("go")
           },
           "da" ={
             dynamics_mainServer("modality")
           }
    )
    
  }, ignoreNULL = FALSE, priority = 1)
  
  observeEvent(input$navBar, {
    currentHash <- sub("#", "", session$clientData$url_hash) # might need to wrap this with `utils::URLdecode` if hash contains encoded characters (not the case here)
    pushQueryString <- paste0("#", input$navBar)
    if(is.null(currentHash) || currentHash != input$navBar){
      freezeReactiveValue(input, "navBar")
      updateQueryString(pushQueryString, mode = "push", session)
    }
  } ,priority=0
  )
  
 #  overview_mainServer("event")
 #  
 #  #Modality analysis
 #  modal_mainServer("modal")
 #  
 # # Differential analysis
 #  
 #  Differential_mainServer("diff")
 #  
 #  #go analysis
 #  go_mainServer("go")
 #  
 #  #dynamics  analysis
 #  
 #  dynamics_mainServer("modality")
  })




























































