#' trackviwer UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
trackviwer_UI<-function(id){
  ns<-NS(id)
  sidebarLayout(
    
    sidebarPanel(
   
      selectInput(ns("TxDb"), label = "Select TxDb package",
                  choices = c("TxDb.Hsapiens.UCSC.hg38.knownGene", 
                              "TxDb.Mmusculus.UCSC.mm10.knownGene")),
      
      selectInput(ns("org"), label = "Select Org package",
                  choices = c("org.Hs.eg.db", 
                              "org.Mm.eg.db")),
      textInput(ns("chr"), label = "chromosome", value = "chr3"),
      numericInput(ns("start"), label = "start", value = 109333900),
      numericInput(ns("end"), label = "end", value = 109338000),
      
      selectInput(
        inputId = ns("theme_choice"),
        label = "Optimize the color theme",
        choices = c("safe", "col", "bw"),  # 你可以扩展这里的主题
        selected = "safe"
      ),
      checkboxInput(ns("trs"), "include transcripts track", value = TRUE),
      
      actionButton(ns("add"), "add coverage track"),
     
      tags$hr(),
      actionButton(ns("refresh"), label="apply change", icon = icon("refresh")),
      actionButton(ns("load"), label="load a saved session", icon = icon("upload")),
      tags$hr(),
      
      
      tags$h4("Set the genomic coordinates for:"),
      actionButton(ns("preSet1"), "Example 1"),
      actionButton(ns("preSet2"), "Example 2"),
      actionButton(ns("preSet3"), "Example 3"),
    
      textInput(ns("symbol"), label="search by gene name", value = ""),
      tags$script(
        'Shiny.addCustomMessageHandler("scrollCallback",
                      function(msg) {
                        window.scrollTo(0, 0);
                      });
                    Shiny.addCustomMessageHandler("loadCallback",
                      function(msg){
                        var id = setInterval(frame, 100);
                        function frame(){
                          if(typeof(d3.select("#importjsonfilename"))!="undefined"){
                             clearInterval(id);
                             document.getElementById("importjsonfilename").click();
                            }
                          }
                        }
                      });
                    ')
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      htmlOutput(ns("geneLabel")),   # ← 新增,
      div(style = "
        padding-left: 80px;     /* 留给标签的空白 */
      ",
          browseTracksOutput(ns("trackViewer"))
      )
    )
  )
}


# Define server logic required to draw tracks
trackviwer_mainServer <- function(id) {
  moduleServer(id, function(input, output,session) {

  global <- reactiveValues(refresh = FALSE, fileIndex=0,
                           fileinserted=c(),obsList=list(), forceRefresh = FALSE)
  observe({
    if(input$refresh){
      isolate(global$refresh <- TRUE)
      isolate(global$forceRefresh <- FALSE)
      session$sendCustomMessage(type="scrollCallback", 1)
      output$trackViewer <- plot
    }else{
      isolate(global$refresh <- FALSE)
      isolate(global$forceRefresh <- FALSE)
    } 
  })
  observeEvent({
    input$TxDb
    input$org
    input$chr
    input$start
    input$end
    input$trs},{
      if(global$forceRefresh){
        isolate(global$refresh <- TRUE)
      }else{
        isolate(global$refresh <- FALSE)
      }
      isolate(global$forceRefresh <- FALSE)
    })
  
  
  
  
  

  
  ## ① 点击  Apply change  时 —— 重新生成标题
  observeEvent(input$refresh, {
    req(input$chr, input$start, input$end, input$TxDb, input$org)
    
    # ❶ 构造区间
    gr <- GRanges(input$chr,
                  IRanges(as.numeric(input$start), as.numeric(input$end)))
    seqlevelsStyle(gr) <- "UCSC"
    
    # ❷ 取 EntrezID
    ids <- tryCatch(
      trackViewer::getGeneIDsFromTxDb(gr, get(input$TxDb)),
      error = function(e) character(0)
    )
    ids <- ids[!is.na(ids)]
    if (!length(ids)) {
      output$geneLabel <- renderUI(NULL);   # 无基因 ⇒ 清空
      return()
    }
    
    # ❸ Entrez ➜ symbol
    symMap <- switch(input$org,
                     "org.Hs.eg.db" = org.Hs.egSYMBOL,
                     "org.Mm.eg.db" = org.Mm.egSYMBOL,
                     org.Hs.egSYMBOL)
    syms <- unique(unlist(mget(ids, symMap, ifnotfound = NA)))
    syms <- syms[!is.na(syms)]
    
    # ❹ 输出蓝色大号标题
    if (length(syms))
      output$geneLabel <- renderUI(
        tags$h3(
          style = "color: #1565c0;          /* 稍深一点的蓝 */
            font-family: 'Helvetica', 'Arial', sans-serif;
            font-weight: 600;        /* 半粗体，比 bold 轻一点 */
            font-size: 20px;         /* 比默认 h3 小一号 */
            line-height: 1.2; margin-top:0;text-align:center;",
          paste(syms, collapse = ", "))
      )
    else
      output$geneLabel <- renderUI(NULL)
  })
  
  ## ② 点击其他任意按钮时 —— 立即清空标题
  observeEvent(
    c(input$add, input$preSet1, input$preSet2, input$preSet3,
      input$TxDb, input$org, input$chr, input$start, input$end),
    {
      output$geneLabel <- renderUI(NULL)    # 覆盖为 NULL ⇒ 页面上消失
    },
    ignoreInit = TRUE       # 不干扰最初界面加载
  )
  
  
  
  
  
  
  
  
  plot <- renderbrowseTracks({
    if(!global$refresh) return()
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message="self checking", value=0)
    gr <- GRanges(input$chr, IRanges(as.numeric(input$start), as.numeric(input$end)))
    gr.NCBI <- gr.UCSC <- gr
    seqlevelsStyle(gr.UCSC) <- "UCSC"
    seqlevelsStyle(gr.NCBI) <- "NCBI"
    progress$set(message="loading library", value=0.03)
    require(input$TxDb, character.only = TRUE)
    require(input$org, character.only = TRUE)
    progress$set(message="get transcripts", value=0.05)
    trs <- tryCatch(geneModelFromTxdb(get(input$TxDb), get(input$org), gr=gr.UCSC), 
                    error = function(e){ NULL })
    progress$set(message="reading track data", value=0.15)
    tks <- list()
    # step = 0.8/(global$fileIndex+global$lolliIndex)
    step= 0.8/global$fileIndex
    if(global$fileIndex>0){
      for(i in seq.int(global$fileIndex)){
        progress$set(message="reading track data", value=0.15+i*step)
        if(paste0("filecontainer", i) %in% global$fileinserted){
          if(paste0("filecontainer", i) %in% global$fileinserted){
            if (input[[paste0("format", i)]] == "bam") {
              files <- input[[paste0("file", i)]]
              if (is.null(files)) {
                showNotification("❌ No BAM/BAI files uploaded.", type = "error")
                return(NULL)
              }
              
              # ⬇️ 识别 BAM 文件和 BAI 文件
              bam_row <- which(grepl("\\.bam$", files$name, ignore.case = TRUE))
              bai_row <- which(grepl("\\.bai$", files$name, ignore.case = TRUE))
              
              if (length(bam_row) != 1) {
                showNotification("❌ Please upload exactly one .bam file.", type = "error")
                return(NULL)
              }
              
              bam_name <- files$name[bam_row]
              bam_path <- files$datapath[bam_row]
              
              expected_bai_path <- paste0(bam_path, ".bai")
              
              # ⬇️ 如果上传了 .bai 文件，重命名为 bam.bai
              if (length(bai_row) == 1) {
                bai_path <- files$datapath[bai_row]
                file.rename(bai_path, expected_bai_path)
              }
              
              # ⬇️ 如果没上传 .bai 文件，就尝试自动生成索引
              if (!file.exists(expected_bai_path)) {
                tryCatch(
                  error = function(e) {
                    showNotification("❌ .bai file missing and indexing failed.", type = "error")
                    return(NULL)
                  })
              }
              
              # ⬇️ 尝试使用 gr.UCSC 读取 BAM 区间
              tks[[input[[paste0("sample", i)]]]] <- tryCatch({
                importBam(
                  file   = bam_path,
                  # pairs  = testPairedEndBam(bam_path),
                  pairs=FALSE,
                  ranges = gr.UCSC
                )
              }, error = function(e1) {
                # ⬇️ 如果 UCSC 样式失败，则尝试 NCBI 样式
                tryCatch({
                  importBam(
                    file   = bam_path,
                    pairs  = testPairedEndBam(bam_path),
                    ranges = gr.NCBI
                  )
                }, error = function(e2) {
                  showNotification("❌ Failed to import BAM file.", type = "error")
                  NULL
                })
              })
              
              
              ## ---------- 读 BAM，先放到 trk ----------
              # trk <- tryCatch({
              #   importBam(
              #     file   = bam_path,
              #     pairs  = FALSE,        # 单端/双端都行，只看 coverage 用 FALSE 最安全
              #     ranges = gr.UCSC
              #   )
              # }, error = function(e1) {
              #   tryCatch({
              #     importBam(
              #       file   = bam_path,
              #       pairs  = FALSE,
              #       ranges = gr.NCBI
              #     )
              #   }, error = function(e2) {
              #     showNotification("❌ Failed to import BAM file.", type = "error")
              #     NULL
              #   })
              # })
              # 
              # ## ---------- 调试输出 ----------
              # if (!is.null(trk)) {
              #   cat("=== BAM import summary ===\n")
              #   cat("Number of ranges :", length(trk@dat), "\n")
              #   cat("Columns in mcols:", names(mcols(trk@dat)), "\n")
              #   
              #   if ("score" %in% names(mcols(trk@dat))) {
              #     cat("First 5 scores  :", head(mcols(trk@dat)$score, 5), "\n")
              #   } else {
              #     cat("⚠️  NO 'score' column found!\n")
              #   }
              #   cat("===========================\n")
              # }
              # 
              # ## ---------- 放回轨道列表 ----------
              # tks[[ input[[paste0("sample", i)]] ]] <- trk
              # 
            }
            
            # if(input[[paste0("format", i)]]=="bam"){
            #   tks[[input[[paste0("sample", i)]]]] <- 
            #     tryCatch(importBam(file = file.path(datafolder, input[[paste0("file", i)]]),
            #                        pairs = testPairedEndBam(file.path(datafolder, input[[paste0("file", i)]])),
            #                        ranges = gr.UCSC),
            #              error = function(e){ 
            #                tryCatch(importBam(file = file.path(datafolder, input[[paste0("file", i)]]),
            #                                   pairs = testPairedEndBam(file.path(datafolder, input[[paste0("file", i)]])),
            #                                   ranges = gr.NCBI),
            #                         error=function(e){NULL})
            #                })
            # }
            
            # 一会修改这个逻辑
            # else{
            #   tks[[input[[paste0("sample", i)]]]] <- 
            #     tryCatch(importScore(file = ifelse(grepl("\\:\\/\\/", input[[paste0("file", i)]]), 
            #                                        input[[paste0("file", i)]],
            #                                        file.path(datafolder, input[[paste0("file", i)]])),
            #                          format = input[[paste0("format", i)]],
            #                          ranges = gr.UCSC), 
            #              error = function(e){ 
            #                tryCatch(importScore(ifelse(grepl("\\:\\/\\/", input[[paste0("file", i)]]), 
            #                                            input[[paste0("file", i)]], 
            #                                            file = file.path(datafolder, input[[paste0("file", i)]])),
            #                                     format = input[[paste0("format", i)]],
            #                                     ranges = gr.NCBI),
            #                         error=function(e){NULL})
            #                })
            # }
          }
        }
      }
    }
    
    trackList <- list()
    if(input$trs && length(trs)>0 && length(tks)==0){
      trackList <- trackList(trs)
    }
    if((length(trs)==0 || !input$trs) && length(tks)>=0){
      trackList <- trackList(tks)
    }
    if((length(trs)>0 && input$trs) && length(tks)>0){
      trackList <- trackList(trs, tks, heightDist = c(1, length(tks)))
    }
    
    progress$set(message="Plot data", value=0.95)
    on.exit(progress$close())
    
    if(length(trackList)>0){
      optSty <- optimizeStyle(trackList, theme=input$theme_choice)
      trackList <- optSty$tracks
      viewerStyle <- optSty$style
      
    
      
      browseTracks(
        trackList,
        gr          = gr.UCSC,
        viewerStyle =viewerStyle
     
        
      )
    }
  })
  output$trackViewer <- plot
  
  
  observeEvent(input$add, {
    isolate(global$fileIndex <- global$fileIndex+1)
    currentIndex <- global$fileIndex
    isolate(global$refresh <- FALSE)
    id = paste0("filecontainer", currentIndex)
    insertUI(
      selector = paste0("#", session$ns("add")),   # ★ 关键改动
      # selector = "#add",
             where = "beforeBegin",
          
             ui <- tags$div(
               tagList(
                 tags$h4("Add data track from BAM file"),
                 
                 # ⬇️ 上传 BAM 和 BAI 文件
                 fileInput(
                   inputId = session$ns(paste0("file", currentIndex)),
                   label   = "Meanwhile upload .bam and .bam.bai files",
                   accept  = c(".bam", ".bai"),
                   multiple = TRUE
                 ),
                 
                 # ⬇️ 格式强制为 bam
                 selectInput(
                   inputId = session$ns(paste0("format", currentIndex)),
                   label   = "file format",
                   choices = c("BAM" = "bam"),
                   selected = "bam"
                 ),
                 
                 # ⬇️ 样本名称输入
                 textInput(
                   inputId = session$ns(paste0("sample", currentIndex)),
                   label   = "sample name",
                   value   = ""
                 ),
                 
                 # ⬇️ 删除按钮
                 actionButton(
                   inputId = session$ns(paste0("remove", currentIndex)),
                   label   = "remove above track",
                   icon    = icon("remove")
                 ),
                 
                 tags$hr()
               ),
               id = id
             )
             
    )
    isolate(global$fileinserted <- c(id, global$fileinserted))
    if(is.null(global$obsList[[paste0(id,"remove")]])){
      isolate(global$obsList[[paste0(id,"remove")]] <- observeEvent(input[[paste0("remove", currentIndex)]],
                                                                    {
                                                                      isolate(global$refresh <- FALSE)
                                                                      removeUI(selector = paste0("#", id))
                                                                      isolate(global$fileinserted <- 
                                                                                global$fileinserted[-which(global$fileinserted==id)])
                                                                    })
      )
    }
    if(is.null(global$obsList[[paste0(id,"fileinput")]])){
      isolate(global$obsList[[paste0(id,"fileinput")]] <- 
                observeEvent(input[[paste0("file", currentIndex)]],
                             {
                               isolate(global$refresh <- FALSE)
                               files <- input[[paste0("file", currentIndex)]]
                               # 1. 获取上传的所有文件信息
                               # files <- input[[paste0("file", i)]]
                               if (is.null(files)) return(NULL)
                               
                               # # 2. 找出 .bam 文件名（长度 = 1）
                               bam_row <- which(grepl("\\.bam$", files$name, ignore.case = TRUE))
                               
                               filename <- files$name[bam_row]
                               filename <- sub(".gz$", "", filename, ignore.case = TRUE)
                               
                               # 3. 使用 tools::file_ext 提取扩展名，并安全匹配
                               fileext_raw <- tolower(tools::file_ext(filename))
                               
                               fileext <- if (length(fileext_raw) == 1 && nzchar(fileext_raw)) {
                                 switch(fileext_raw,
                                        "bam"      = "bam",
                                        "unknown")
                               }
                               
                         
                               if(fileext!="unknown"){
                                 updateSelectInput(session, inputId = paste0("format", currentIndex),
                                                   selected =fileext )
                               }
                               filename <- sub(paste0(".", fileext, "$"), "", filename)
                               updateTextInput(session, inputId = paste0("sample", currentIndex),
                                               value = filename)
                             })
      )
    }
  })
  
  
  observeEvent(input$load, {
    output$trackViewer <- 
      renderbrowseTracks({
        A=new("track", dat=GRanges(1, IRanges(c(1, 3), c(2, 4)), score=c(1, 10)), 
              type="data", format="BED")
        browseTracks(trackList(A), gr=GRanges(1, IRanges(1, 10000)))
      })
    session$sendCustomMessage(type="loadCallback", 1)
  })
  
  observeEvent(input$preSet1, {
    updateSelectInput(session, inputId = "TxDb", selected = "TxDb.Hsapiens.UCSC.hg38.knownGene")
    updateSelectInput(session, inputId = "org", selected = "org.Hs.eg.db")
    updateTextInput(session, inputId = "chr", value = "chr1")
    updateNumericInput(session, inputId = "start", value = 182839300)
    updateNumericInput(session, inputId = "end", value = 182842862)
  })
  observeEvent(input$preSet2, {
    updateSelectInput(session, inputId = "TxDb", selected = "TxDb.Hsapiens.UCSC.hg38.knownGene")
    updateSelectInput(session, inputId = "org", selected = "org.Hs.eg.db")
    updateTextInput(session, inputId = "chr", value = "chr3")
    updateNumericInput(session, inputId = "start", value = 109333900)
    updateNumericInput(session, inputId = "end", value = 109338000)
  })
  observeEvent(input$preSet3, {
    updateSelectInput(session, inputId = "TxDb", selected = "TxDb.Hsapiens.UCSC.hg38.knownGene")
    updateSelectInput(session, inputId = "org", selected = "org.Hs.eg.db")
    updateTextInput(session, inputId = "chr", value = "chr6")
    updateNumericInput(session, inputId = "start", value = 136277843)
    updateNumericInput(session, inputId = "end", value = 136279892)
  })

  observeEvent(input$symbol, {
    isolate(global$refresh <- FALSE)
    if(input$symbol!=""){
      progress <- shiny::Progress$new()
      progress$set(message=paste("searching", input$symbol), value=0)
      require(input$TxDb, character.only = TRUE)
      require(input$org, character.only = TRUE)
      prefix <- sub(".db", "", input$org)
      eid <- tryCatch(mget(input$symbol, envir=get(paste0(prefix, "ALIAS2EG")), ifnotfound = NA),
                      error=function(e){
                        message(e)
                        return(NA)
                      })
      eid <- eid[[1]]
      progress$set(message=paste("get gene id", eid), value=50)
      if(!is.na(eid)){
        location <- tryCatch(AnnotationDbi::select(get(input$TxDb), keys=eid, 
                                    columns=c("EXONCHROM", "EXONSTART", "EXONEND"), 
                                    keytype = "GENEID"), 
                             error=function(e){
                               message(e)
                               return(data.frame())
                             })
        progress$set(message="get gene location", value=90)
        if(nrow(location)>0){
          seqnames <- as.character(location$EXONCHROM)[1]
          rg <- range(c(location$EXONSTART, location$EXONEND))
          progress$set(message="set genomic loation", value=0.95)
          updateTextInput(session, inputId = "chr", value = seqnames)
          updateNumericInput(session, inputId = "start", value = rg[1])
          updateNumericInput(session, inputId = "end", value = rg[2])
        }
      }
      session$sendCustomMessage(type="scrollCallback", 1)
      on.exit(progress$close())
    }
  })
  }
  )
}
