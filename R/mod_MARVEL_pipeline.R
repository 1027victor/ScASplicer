#' MARVEL_pipline UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_MARVEL_pipeline_ui <- function(id){
  ns <- NS(id)
  tagList(
  fluidRow(
    column(2,
           wellPanel(
             h4("Settings"),
             numericInput(ns("coverage_threshold"), p("Set Coverage Threshold",style="color:black"), value = 10, min = 1, step = 1)
             
           )
    ),
    
    column(6,
           wellPanel(
             h4("File Inputs"),
             fluidRow(
               column(6, fileInput(ns("file_pheno"), "Sample metadata", accept = c(".txt"))),
               column(6, fileInput(ns("file_sj"), "Splice junction counts matrix", accept = c(".txt")))
             ),
             fluidRow(
               column(6, fileInput(ns("file_se"), "SE_featureData", accept = c(".txt"))),
               column(6, fileInput(ns("file_mxe"), "MXE_featureData", accept = c(".txt")))
             ),
             fluidRow(
               column(6, fileInput(ns("file_ri"), "RI_featureData", accept = c(".txt"))),
               column(6, fileInput(ns("file_a5ss"), "A5SS_featureData", accept = c(".txt")))
             ),
             fluidRow(
               column(6, fileInput(ns("file_a3ss"), "A3SS_featureData", accept = c(".txt"))),
               column(6, fileInput(ns("file_intron"), "Intron count matrix", accept = c(".txt")))
             ),
             fluidRow(
               column(6, fileInput(ns("file_tpm"), "Gene expression matrix", accept = c(".txt"))),
               column(6, fileInput(ns("file_tpm_feature"), "Gene metadata", accept = c(".txt")))
             ),
             fluidRow(
               column(12, fileInput(ns("file_gtf"), "Gene transfer file (GTF)"))
             )
           ),offset = 1
    )
    
  ),
  br(),
  fluidRow(
    column(3),
    column(3,actionButton(ns("run_pipeline"), "Run Pipeline", class='btn-info',icon = icon("running"))),
    column(3,downloadButton(ns("downloadData"), class='btn-primary',"Download MARVEL Object")),
    column(3))
  )
}
    
#' MARVEL_pipline Server Functions
#'
#' @noRd 
mod_MARVEL_pipeline_server <- function(id){
  moduleServer( id, function(input, output, session){
    marvel_data <- reactiveVal(NULL)
    
    observeEvent(input$run_pipeline, {

      if(is.null(input$file_pheno) || is.null(input$file_sj) || is.null(input$file_se) || is.null(input$file_mxe) || is.null(input$file_ri) || is.null(input$file_a5ss) || is.null(input$file_a3ss) || is.null(input$file_intron) || is.null(input$file_tpm) || is.null(input$file_tpm_feature) || is.null(input$file_gtf)){
        showModal(modalDialog(
          title = "Error",
          "Please upload all files",
          footer = modalButton("close")
        ))
        return(NULL)
      }
      
      
      # 读取上传的文件
    withProgress(message = 'processing...', value = 0, {
      df.pheno <- read.table(input$file_pheno$datapath, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA")
      sj <- as.data.frame(fread(input$file_sj$datapath, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA"))
      df.feature.se <- read.table(input$file_se$datapath, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA")
      df.feature.mxe <- read.table(input$file_mxe$datapath, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA")
      df.feature.ri <- read.table(input$file_ri$datapath, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA")
      df.feature.a5ss <- read.table(input$file_a5ss$datapath, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA")
      df.feature.a3ss <- read.table(input$file_a3ss$datapath, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA")
      
      df.feature.list <- list(df.feature.se, df.feature.mxe, df.feature.ri, df.feature.a5ss, df.feature.a3ss)
      names(df.feature.list) <- c("SE", "MXE", "RI", "A5SS", "A3SS")
      
      df.intron.counts <- as.data.frame(fread(input$file_intron$datapath, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA"))
      df.tpm <- read.table(input$file_tpm$datapath, sep="\t", header=TRUE, stringsAsFactors=FALSE)
      df.tpm.feature <- read.table(input$file_tpm_feature$datapath, sep="\t", header=TRUE, stringsAsFactors=FALSE)
      gtf <- as.data.frame(fread(input$file_gtf$datapath, sep="\t", header=FALSE, stringsAsFactors=FALSE, na.strings="NA", quote="\""))
      
      # 创建 Marvel 对象
      marvel <- CreateMarvelObject(SpliceJunction=sj,
                                   SplicePheno=df.pheno,
                                   SpliceFeature=df.feature.list,
                                   IntronCounts=df.intron.counts,
                                   GeneFeature=df.tpm.feature,
                                   Exp=df.tpm,
                                   GTF=gtf
      )
      
      # 使用用户设置的 Coverage Threshold 进行 PSI 计算
      coverage_threshold <- input$coverage_threshold
      marvel <- ComputePSI(MarvelObject=marvel, CoverageThreshold=coverage_threshold, UnevenCoverageMultiplier=10, EventType="SE")
      marvel <- ComputePSI(MarvelObject=marvel, CoverageThreshold=coverage_threshold, UnevenCoverageMultiplier=10, EventType="MXE")
      marvel <- ComputePSI(MarvelObject=marvel, CoverageThreshold=coverage_threshold, EventType="RI", thread=4)
      marvel <- ComputePSI(MarvelObject=marvel, CoverageThreshold=coverage_threshold, EventType="A5SS")
      marvel <- ComputePSI(MarvelObject=marvel, CoverageThreshold=coverage_threshold, EventType="A3SS")
      
      
      index.1 <- which(df.pheno$cell.type %in% unique(df.pheno$cell.type))
      index.2 <- which(df.pheno$qc.seq=="pass")
      
      # sample IDs
      index <- Reduce(intersect, list(index.1, index.2))
      sample.ids <- df.pheno[index, "sample.id"]
      
      marvel <- SubsetSamples(MarvelObject=marvel,
                              sample.ids=sample.ids
      )
      
      marvel <- TransformExpValues(MarvelObject=marvel,
                                   offset=1,
                                   transformation="log2",
                                   threshold.lower=1
      )
      
      marvel <- CheckAlignment(MarvelObject=marvel, level="splicing")
      
      # Check gene data
      marvel <- CheckAlignment(MarvelObject=marvel, level="gene")
      
      # Cross-check splicing and gene data
      marvel <- CheckAlignment(MarvelObject=marvel, level="splicing and gene")
      
      showModal(modalDialog(
        title = "",
        "The data processing is complete!",
        footer = modalButton("close")
      ))
      
      marvel_data(marvel)
      

    }
    )
      })
    
    # 下载按钮
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("result", Sys.Date(), ".RData", sep="")
      },
      content = function(file) {
        marvel <- marvel_data()
        if (!is.null(marvel)) {
          save(marvel, file=file)
        }
      }
    )
 
  })
}
    
## To be copied in the UI
# mod_MARVEL_pipline_ui("MARVEL_pipline_1")
    
## To be copied in the server
# mod_MARVEL_pipline_server("MARVEL_pipline_1")
