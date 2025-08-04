#' The application server-side
#'
#' @param input,output Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import MARVEL
#' @import rtracklayer
#' @import VariantAnnotation
#' @import Rsamtools
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import TxDb.Mmusculus.UCSC.mm10.knownGene
#' @import GenomicRanges
#' @import IRanges
#' @import ggnewscale
#' @import ggrepel
#' @import parallel
#' @import reshape2
#' @import stringr
#' @import textclean
#' @import ggplot2
#' @import gridExtra
#' @import highcharter
#' @import shinycssloaders
#' @import dplyr
#' @import scales
#' @import tidyr
#' @import circlize
#' @import grid
#' @import graphics
#' @import ComplexHeatmap
#' @import factoextra
#' @import FactoMineR
#' @import fitdistrplus
#' @import kSamples
#' @import twosamples
#' @import AnnotationDbi
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import Biostrings
#' @import BSgenome
#' @import BSgenome.Hsapiens.NCBI.GRCh38
#' @import data.table
#' @noRd

# source function
source("R/utils.R")

app_server<-function(input, output,session) {
  options(shiny.maxRequestSize=300*1024^3)
  stopifnot(packageVersion("trackViewer")>="1.19.11")
  pipeline_server <- mod_MARVEL_pipeline_server("pipeline")
  event_server <- overview_mainServer("event")
  modal_server <- modal_mainServer("modal")
  diff_server <- Differential_mainServer("diff")
  go_server <- go_mainServer("go")
  modality_server <- dynamics_mainServer("modality")
  vract_server<-trackviwer_mainServer("vrcat")
  
  observeEvent(session$clientData$url_hash, {
    currentHash <- sub("#", "", session$clientData$url_hash)
    if(is.null(input$navBar) || !is.null(currentHash) && currentHash != input$navBar){
      freezeReactiveValue(input, "navBar")
      updateNavbarPage(session, "navBar", selected = currentHash)
    }
    
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
}





























































