#' goUsage_main UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
# go_ui
go_UI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      
      radioButtons(ns("selection"),'',c(input_Differential_gene_analysis_results='input Differential gene analysis results',input_target_gene_file='input target gene file'),selected = 'input Differential gene analysis results',inline=T),
      uiOutput(ns("selection_ui"))
      
      
    ),
    
    mainPanel(splitLayout(div(helpText("Show the results in this view.",style="font-weight:bold;font-size:24px;color:#303133 !important;font-family: OpenSans-Semibold, Helvetica Neue, Helvetica, Arial, sans-serif !important;"),
                              uiOutput(ns("go_plot_scroll_box")),helpText("Show output data.",style="font-weight:bold;font-size:24px;color:#303133 !important;font-family: OpenSans-Semibold, Helvetica Neue, Helvetica, Arial, sans-serif !important;"),br(),uiOutput(ns('go_download')),br(),DT::dataTableOutput(ns('go_result'))),  wellPanel(div(numericInput(ns('go_term'),'GO_term',10,1,20),uiOutput(ns("plot_setting")))),
                          cellWidths = c("70%","30%"))
    )
  )
}

#go_server
go_mainServer <- function(id) {
  moduleServer(id, function(input, output,session) {
    observe({
      if (input$selection=='input Differential gene analysis results')
      {
        output$selection_ui<-renderUI({
          div(fileInput(session$ns("go_file"), label =p("File input:",style="color:black; text-align:center")),
              radioButtons(session$ns("sep"),'Sep',c(Tab='\t',Comma=',',Semicolon=';'),selected = '\t',inline=T),
              sliderTextInput(session$ns("pvalue_slider"),
                              label ='p.val.adj',
                              grid = TRUE,
                              choices = c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.3, 0.5, 1),
                              selected = 0.1),
              fluidRow(column(6,textInput(session$ns("log2fc"),"log2fc",value = 7)),
                       column(6,textInput(session$ns("neglog2fc"),"-log2fc",value = -5))),
              radioButtons(session$ns("species"),'species',c(human='org.Hs.eg.db',mouse='org.Mm.eg.db'),selected = 'org.Hs.eg.db',inline=T),             
              radioButtons(session$ns("caterogy"),'Caterogy',c(BP='BP',CC='CC',MF='MF'),selected = 'BP',inline=T),
              actionButton(session$ns("Gene_ontology_analysis"), "go analysis",class = "btn-info"),
              selectInput(session$ns("tools"),"Choose Bubble char or Bar chart or Circle diagram?",choices = c("Bubble char","Bar chart","Circle diagram")),
              actionButton(session$ns("go_action"), "Start Plot",class = "btn-primary")
          )
        }
        )
      }
      
      else if(input$selection=='input target gene file')
      {
        output$selection_ui<-renderUI({
          div(fileInput(session$ns("goal_go_file"), label =p("File input:",style="color:black; text-align:center")),
              radioButtons(session$ns("goal_species"),'species',c(human='org.Hs.eg.db',mouse='org.Mm.eg.db'),selected = 'org.Hs.eg.db',inline=T),
              radioButtons(session$ns("goal_caterogy"),'Caterogy',c(BP='BP',CC='CC',MF='MF'),selected = 'BP',inline=T),
              sliderTextInput(session$ns("pvalueCutoff"),
                              label ='pvalueCutoff',
                              grid = TRUE,
                              choices = c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.3, 0.5, 1),
                              selected = 0.05),
              sliderTextInput(session$ns("qvalueCutoff"),
                              label ='qvalueCutoff',
                              grid = TRUE,
                              choices = c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.3, 0.5, 1),
                              selected = 0.05),    
              fluidRow(column(6,textInput(session$ns("minGSSize"),"minGSSize",value = 10)),
                       column(6,textInput(session$ns(" maxGSSize")," maxGSSize",value = 500))),
              actionButton(session$ns("goal_Gene_ontology_analysis"), "go analysis",class = "btn-info"),
              selectInput(session$ns("goal_tools"),"Choose Bubble char or Bar chart or Circle diagram?",choices = c("Bubble char","Bar chart")),
              actionButton(session$ns("goal_go_action"), "Start Plot",class = "btn-primary")
          )
        }
        )
      }
    }
    )
    
    
    
    
    observeEvent(input$tools, {
      if (input$tools == "Circle diagram") {
        output$plot_setting<-renderUI({
          div(colourpicker::colourInput(session$ns("col1_circle"),"Outer color","#69c3c5"),
              colourpicker::colourInput(session$ns("col2_circle"),"Outer font color","#FFFFFF"),
              colourpicker::colourInput(session$ns("col3_circle"),"Upregulated","#a0c6c9"),
              colourpicker::colourInput(session$ns("col4_circle"),"Downregulated","#c05678"),
              colourpicker::colourInput(session$ns("col5_circle"),"Fourth circle color","#69c3c5"),
              fluidRow(column(6,colourpicker::colourInput(session$ns("col6_circle"),"log10_legend","#FF906F")),
                       column(6,colourpicker::colourInput(session$ns("col7_circle"),"","#861D30"))),     
              fluidRow(column(6,textInput(session$ns("plotcirclewidth"),"Plot width",value = 10)),
                       column(6,textInput(session$ns("plotcircleheight"),"Plot height",value = 10))),
              radioButtons(session$ns("extPlot"), 'Plot output format',choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
              downloadButton(session$ns("plotdownloadData"),"Download Plot"))
          
        })
      }
      else if(input$tools == "Bar chart"){
        output$plot_setting<-renderUI({
          div(colourpicker::colourInput(session$ns("col1_bar"),"color min","#e1e8f0"),
              colourpicker::colourInput(session$ns("col2_bar"),"color max","#326fa4"),
              textInput(session$ns("fontSize"),"Font size","12px"),
              fluidRow(
                column(6,textInput(session$ns("plotbarWidth"), "Plot Width", value = 850)),
                column(6,textInput(session$ns("plotbarHeight"), "Plot Height",value = 600))
              ))
          
        })
      }
      
      else if(input$tools == "Bubble char"){
        output$plot_setting<-renderUI({
          div(colourpicker::colourInput(session$ns("col1_bb"),"color min","#FFF6A5"),
              colourpicker::colourInput(session$ns("col2_bb"),"color max","#c1021f"),
              textInput(session$ns("Bubble_fontSize"),"Font size","12px"),
              fluidRow(
                column(6,textInput(session$ns("mixsize"), "Min size", value = 10)),
                column(6,textInput(session$ns("maxsize"), "Max size", value = 30),
                )),
              fluidRow(
                column(6,textInput(session$ns("plotbubbleWidth"), "Plot Width", value = 850)),
                column(6,textInput(session$ns("plotbubbleHeight"), "Plot Height",value = 600))
              )
          )
          
        })
      }
    })
    
    
    go_data<-reactiveValues()
    #load data
    observeEvent(input$Gene_ontology_analysis,{
      
      withProgress(message = 'processing...', value = 0, {
        go_data$go=NULL
        go_data$circle_go=NULL
        infile <- input$go_file
        req(infile)
        # if (is.null(infile)) {
        #     return(NULL)
        # }
        # else {
        
        log2fc <- as.numeric(input$log2fc)
        neglog2fc <- as.numeric(input$neglog2fc)
        p.val.adj <- as.numeric(input$pvalue_slider)
        # withProgress(message = 'processing...', value = 0, {
        data <- read.csv(infile$datapath,sep=input$sep, header = TRUE)
        data<-data[data$log2fc>log2fc|data$log2fc<neglog2fc,]
        data<-data[data$p.val.adj<p.val.adj,]
        # 找出第二列中非重复（即第一次出现的）行
        non_duplicated_rows <- !duplicated(data[, 2])
        
        # 保留这些非重复行
        data_unique <- data[non_duplicated_rows, ]
        rownames(data_unique) <- as.character(data_unique[, 2])
        data<-data_unique[,-2]
        up <- rownames(data)[data$sig == "up"]#差异上调
        down <- rownames(data)[data$sig == "down"]#差异下调
        diff <- c(up, down)#所有差异基因
        diff_entrez <- bitr(diff,
                            fromType = "SYMBOL",
                            toType = "ENTREZID",
                            OrgDb = input$species)
        # head(diff_entrez)
        
        #单个Ontology富集(CC为例):
        GO_diff <- enrichGO(gene = diff_entrez$ENTREZID,
                            OrgDb = input$species, #指定包含该物种注释信息的org包
                            ont = input$caterogy, #"BP"，"MF","CC"三选一,或"ALL"合并
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            minGSSize = 10,
                            maxGSSize = 500,
                            readable = T) #结果表中将entrez转回symbol
        # browser()
        #提取结果表格：
        GO_result <- GO_diff@result
        GO_result <- GO_diff@result %>%
          mutate(neg_log10pvalue = -log10(pvalue))%>%
          arrange(desc(neg_log10pvalue)) 
        # head(10)
        # head(input$go_term)
        # GO_result$color <- scales::col_numeric(palette = c(input$col1bb, input$col2bb), 
        #                              domain = range(GO_result$neg_log_pvalue))( GO_result$neg_log_pvalue)
        #delete infinite rows
        infinite_rows <- is.infinite(GO_result$neg_log10pvalue)
        GO_result <- GO_result[!infinite_rows, ]
        
        data_one<-data[,'sig', drop = FALSE]
        data_one['gene_id']<-row.names(data)
        row.names(data_one)<-NULL
        GO_result_split <- GO_result %>%separate_rows(geneID, sep = "/")
        data_one <- rename(data_one, geneID=gene_id)
        merged_data <- merge(GO_result_split, data_one, by = "geneID")   
        proportions <- merged_data %>%
          group_by(ID) %>%
          summarise(
            Total = n(),
            Upregulated = sum(sig == "up"),
            Downregulated = sum(sig == "down"),
            Upregulated_Prop = Upregulated / Total,
            Downregulated_Prop = Downregulated / Total
          )
        merge_all<-merge(GO_result,proportions,by = "ID")
        rownames(merge_all)<-merge_all$ID
        merge_all$gene_num.min <- 0
        merge_all$gene_num.max <- as.numeric(strsplit(merge_all$BgRatio[1], "/")[[1]][2])
        merge_all$gene_num.rich <- as.numeric(unlist(lapply(merge_all$BgRatio, 
                                                            function(x) strsplit(x,"/")[[1]][1])))
        rich_gene_num <- as.numeric(unlist(lapply(merge_all$GeneRatio,
                                                  function(x) strsplit(x,"/")[[1]][1])))
        # merge_all$"-log10Pvalue" <- -log10(merge_all$pvalue)
        merge_all$rich.factor <- rich_gene_num/merge_all$gene_num.rich
        merge_all$category <-input$caterogy
        merge_all <- merge_all[order(merge_all$p.adjust), ]
        
        setProgress(message = sprintf("Processing progress"))
        # showNotification("complete! ", type = "warning") 
        
        showModal(modalDialog(
          title = "",
          "The data processing is complete!",
          footer = modalButton("close")
        ))
        
        go_data$go=GO_result
        go_data$circle_go=merge_all
        
      })
      # }
    },ignoreInit = TRUE)
    
    
    
    values <- reactiveValues()
    
    
    #start to plot and download
    observeEvent(input$go_action, {
      
      values$plotbarWidth <- input$plotbarWidth
      values$plotbarHeight <- input$plotbarHeight
      values$plotbubbleWidth<-input$plotbubbleWidth
      values$plotbubbleHeight<-input$plotbubbleHeight
      values$plotcirclewidth<-input$plotcirclewidth
      values$plotcircleheight<-input$plotcircleheight
      values$go_term <- input$go_term
      values$col1_bar <- input$col1_bar
      values$col2_bar <- input$col2_bar
      values$col1_bb <- input$col1_bb
      values$col2_bb <- input$col2_bb
      values$Bubble_fontSize <- input$Bubble_fontSize
      values$mixsize <- input$mixsize
      values$maxsize <- input$maxsize
      values$fontSize <- input$fontSize
      values$col1_circle <- input$col1_circle
      values$col2_circle <- input$col2_circle
      values$col3_circle <- input$col3_circle
      values$col4_circle <- input$col4_circle
      values$col5_circle <- input$col5_circle
      values$col6_circle <- input$col6_circle
      values$col7_circle <- input$col7_circle
      # browser()
      if (input$tools=="Bar chart")
      {
        output$go_plot_scroll_box <- renderUI({
          if (is.null(go_data$go)){
            return (NULL)
          }
          else
          {div(style ="overflow-y: scroll; overflow-x: scroll;max-height:700px;border: 1px solid #ddd; padding: 20px;",
               highchartOutput(session$ns("go_plot"),width=values$plotbarWidth,height=values$plotbarHeight)
               #  highchartOutput(session$ns("go_plot"))
               # highchartOutput("go_plot",width="1500px",height="800px")
          )}
        })
        # browser()
        output$go_plot <- renderHighchart({
          my_bar(go_data$go,go_term = values$go_term ,color_min = values$col1_bar, color_max = values$col2_bar,
                 fontSize=values$fontSize,sourceWidth=values$plotbarWidth,sourceHeight=values$plotbarHeight)
          # my_bar(go_data())
        })
        
        output$go_download <- renderUI({
          if (is.null(go_data$go)) {
            return(NULL)
          }
          div(class = "dropdown",
              tags$button(class = "btn btn-custom dropdown-toggle", type = "button", id = "dropdownMenuButton",
                          `data-bs-toggle` = "dropdown", `aria-haspopup` = "true", `aria-expanded` = "false",
                          "Download  Data ", span(class = "caret")
              ),
              div(class = "dropdown-menu", `aria-labelledby` = "dropdownMenuButton",
                  downloadLink(session$ns("download_go_csv"), "CSV", class = "btn-link"),
                  downloadLink(session$ns("download_go_txt"), "TXT", class = "btn-link")
              ))
          
        })
        
        
        output$go_result <- DT::renderDataTable({
          if (is.null(go_data$go)){
            return(NULL)}
          DT::datatable(go_data$go,
                        options = list(pageLength =10, scrollX = TRUE))}
        )
        
        output$download_go_csv <- downloadHandler(
          filename = function() {
            paste("result", Sys.Date(), ".csv", sep="")
          },
          content = function(file) {
            # load your_data_frame 
            write.csv(go_data$go, file)
          }
        )
        
        # downlaoad txt
        output$download_go_txt <- downloadHandler(
          filename = function() {
            paste("result", Sys.Date(), ".txt", sep="")
          },
          content = function(file) {
            # load your_data_frame
            write.table(go_data$go, file, row.names = FALSE)
          }
        )
        
      }
      else if (input$tools=="Bubble char") {
        output$go_plot_scroll_box <- renderUI({
          if (is.null(go_data$go)){
            return (NULL)
          }
          else
          {div(style ="overflow-y: scroll; overflow-x: scroll;max-height:700px;border: 1px solid #ddd; padding: 20px;",
               highchartOutput(session$ns("go_plot"),width=values$plotbubbleWidth,height=values$plotbubbleHeight)
               #  highchartOutput(session$ns("go_plot"))
               # highchartOutput("go_plot",width="1500px",height="800px")
          )}
        })
        
        output$go_plot <- renderHighchart({
          my_Bubble_plot(go_data$go,go_term = values$go_term,color_min = values$col1_bb, color_max = values$col2_bb,
                         sourceWidth=values$plotbubbleWidth,sourceHeight=values$plotbubbleHeight,
                         fontSize=values$Bubble_fontSize,minsize=values$mixsize,maxsize=values$maxsize)
        })
        
        output$go_download <- renderUI({
          if (is.null(go_data$go)) {
            return(NULL)
          }
          div(class = "dropdown",
              tags$button(class = "btn btn-custom dropdown-toggle", type = "button", id = "dropdownMenuButton",
                          `data-bs-toggle` = "dropdown", `aria-haspopup` = "true", `aria-expanded` = "false",
                          "Download  Data ", span(class = "caret")
              ),
              div(class = "dropdown-menu", `aria-labelledby` = "dropdownMenuButton",
                  downloadLink(session$ns("download_go_bubble_csv"), "CSV", class = "btn-link"),
                  downloadLink(session$ns("download_go_bubble_txt"), "TXT", class = "btn-link")
              ))
          
        })
        
        
        output$go_result <- DT::renderDataTable({
          if (is.null(go_data$go)){
            return(NULL)}
          DT::datatable(go_data$go,
                        options = list(pageLength =10, scrollX = TRUE))}
        )
        
        output$download_go_bubble_csv <- downloadHandler(
          filename = function() {
            paste("result", Sys.Date(), ".csv", sep="")
          },
          content = function(file) {
            # load your_data_frame 
            write.csv(go_data$go, file)
          }
        )
        
        # downlaoad txt
        output$download_go_bubble_txt <- downloadHandler(
          filename = function() {
            paste("result", Sys.Date(), ".txt", sep="")
          },
          content = function(file) {
            # load your_data_frame
            write.table(go_data$go, file, row.names = FALSE)
          }
        )
        
      }
      else if (input$tools=="Circle diagram"){
        output$go_plot_scroll_box <- renderUI({
          if (is.null(go_data$circle_go)){
            return (NULL)
          }
          else
          {div(style ="overflow-y: scroll; overflow-x: scroll;max-height:1000px;border: 1px solid #ddd; padding: 20px;",
               plotOutput(session$ns("plot"),width="800px",height="800px")
               # highchartOutput("go_plot",width="1500px",height="800px")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
          )}
        })
        
        output$plot <- renderPlot({
          my_circle(go_data$circle_go, go_term = values$go_term,
                    outcolor = values$col1_circle, outfontcolor = values$col2_circle,
                    upregulated = values$col3_circle, downregulated = values$col4_circle,
                    Fourthcirclecolor = values$col5_circle,lengend_down = values$col6_circle,
                    lengend_up = values$col7_circle)
        })
        
        
        #download
        
        output$go_download <- renderUI({
          if (is.null(go_data$circle_go)) {
            return(NULL)
          }
          div(class = "dropdown",
              tags$button(class = "btn btn-custom dropdown-toggle", type = "button", id = "dropdownMenuButton",
                          `data-bs-toggle` = "dropdown", `aria-haspopup` = "true", `aria-expanded` = "false",
                          "Download  Data ", span(class = "caret")
              ),
              div(class = "dropdown-menu", `aria-labelledby` = "dropdownMenuButton",
                  downloadLink(session$ns("download_circle_go_csv"), "CSV", class = "btn-link"),
                  downloadLink(session$ns("download_circle_go_txt"), "TXT", class = "btn-link")
              ))
          
        })
        
        
        output$go_result <- DT::renderDataTable({
          if (is.null(go_data$circle_go)){
            return(NULL)}
          DT::datatable(go_data$circle_go,
                        options = list(pageLength =10, scrollX = TRUE))}
        )
        
        output$download_circle_go_csv <- downloadHandler(
          filename = function() {
            paste("result", Sys.Date(), ".csv", sep="")
          },
          content = function(file) {
            # load your_data_frame 
            write.csv(go_data$circle_go, file)
          }
        )
        
        # downlaoad txt
        output$download_circle_go_txt <- downloadHandler(
          filename = function() {
            paste("result", Sys.Date(), ".txt", sep="")
          },
          content = function(file) {
            # load your_data_frame
            write.table(go_data$circle_go, file, row.names = FALSE)
          }
        )
        
        # download plot
        
        output$plotdownloadData<-downloadHandler(
          filename = function() {
            paste('circleplot', Sys.Date(), '.',input$extPlot, sep='')
          },
          content = function(file) {
            if (input$extPlot == "pdf") {
              pdf(file,width=as.numeric(values$plotcirclewidth),height=as.numeric(values$plotcircleheight))
            } else if (input$extPlot == "jpeg") {
              jpeg(file, width = as.numeric(values$plotcirclewidth), height = as.numeric(values$plotcircleheight),units = "in",res = 300)
            } else if (input$extPlot == "png") {
              png(file, width = as.numeric(values$plotcirclewidth), height = as.numeric(values$plotcircleheight),units = "in",res = 300)
            }
            
            # plot
            my_circle(go_data$circle_go, go_term = values$go_term,
                      outcolor = values$col1_circle, outfontcolor = values$col2_circle,
                      upregulated = values$col3_circle, downregulated = values$col4_circle,
                      Fourthcirclecolor = values$col5_circle,lengend_down = values$col6_circle,
                      lengend_up = values$col7_circle)
            
            # close
            dev.off()
          }
        )
        
        
        
      }
      
      
      
      
      
      
    }
    )
    
    observeEvent(input$goal_tools, {
      
      if(input$goal_tools == "Bar chart"){
        output$plot_setting<-renderUI({
          div(colourpicker::colourInput(session$ns("goal_col1_bar"),"color min","#e1e8f0"),
              colourpicker::colourInput(session$ns("goal_col2_bar"),"color max","#326fa4"),
              textInput(session$ns("goal_fontSize"),"Font size","12px"),
              fluidRow(
                column(6,textInput(session$ns("goal_plotbarWidth"), "Plot Width", value = 850)),
                column(6,textInput(session$ns("goal_plotbarHeight"), "Plot Height",value = 600))
              ))
          
        })
      }
      
      else if(input$goal_tools == "Bubble char"){
        output$plot_setting<-renderUI({
          div(colourpicker::colourInput(session$ns("goal_col1_bb"),"color min","#FFF6A5"),
              colourpicker::colourInput(session$ns("goal_col2_bb"),"color max","#c1021f"),
              textInput(session$ns("goal_Bubble_fontSize"),"Font size","12px"),
              fluidRow(
                column(6,textInput(session$ns("goal_mixsize"), "Min size", value = 10)),
                column(6,textInput(session$ns("goal_maxsize"), "Max size", value = 30),
                )),
              fluidRow(
                column(6,textInput(session$ns("goal_plotbubbleWidth"), "Plot Width", value = 850)),
                column(6,textInput(session$ns("goal_plotbubbleHeight"), "Plot Height",value = 600))
              )
          )
          
        })
      }
    }
    )
    
    goal<-reactiveValues()
    #load data
    observeEvent(input$goal_Gene_ontology_analysis,{
      withProgress(message = 'processing...', value = 0, {
        goal$goal_go=NULL
        goal_infile <- input$goal_go_file
        if (is.null(goal_infile)) {
          return(NULL)
        }
        # withProgress(message = 'processing...', value = 0, {
        data <- read.csv(goal_infile$datapath,header = TRUE)
        # data=read.xlsx(goal_infile$datapath,startRow = 1)
        diff_entrez <- bitr(data[,],
                            fromType = "SYMBOL",
                            toType = "ENTREZID",
                            OrgDb = input$goal_species)
        # head(diff_entrez)
        
        #单个Ontology富集(CC为例):
        GO_diff <- enrichGO(gene = diff_entrez$ENTREZID,
                            OrgDb = input$goal_species, #指定包含该物种注释信息的org包
                            ont = input$goal_caterogy, #"BP"，"MF","CC"三选一,或"ALL"合并
                            pAdjustMethod = "BH",
                            pvalueCutoff = input$pvalueCutoff,
                            qvalueCutoff = input$qvalueCutoff,
                            minGSSize = input$minGSSize,
                            maxGSSize = input$maxGSSize,
                            readable = T) #结果表中将entrez转回symbol
        # browser()
        #提取结果表格：
        GO_result <- GO_diff@result
        GO_result <- GO_diff@result %>%
          mutate(neg_log10pvalue = -log10(pvalue))%>%
          arrange(desc(neg_log10pvalue)) 
        
        
        
        #delete infinite rows
        infinite_rows <- is.infinite(GO_result$neg_log10pvalue)
        GO_result <- GO_result[!infinite_rows, ]
        
        setProgress(message = sprintf("Processing progress"))
        showNotification("complete! ", type = "warning") 
        
        goal$goal_go=GO_result
        
      })
    },ignoreInit = TRUE)
    
    
    
    goal_values <- reactiveValues()
    
    
    #start to plot and download
    observeEvent(input$goal_go_action, {
      
      goal_values$goal_plotbarWidth <- input$goal_plotbarWidth
      goal_values$goal_plotbarHeight <- input$goal_plotbarHeight
      goal_values$goal_plotbubbleWidth<-input$goal_plotbubbleWidth
      goal_values$goal_plotbubbleHeight<-input$goal_plotbubbleHeight
      goal_values$go_term <- input$go_term
      goal_values$goal_col1_bar <- input$goal_col1_bar
      goal_values$goal_col2_bar <- input$goal_col2_bar
      goal_values$goal_col1_bb <- input$goal_col1_bb
      goal_values$goal_col2_bb <- input$goal_col2_bb
      goal_values$goal_Bubble_fontSize <- input$goal_Bubble_fontSize
      goal_values$goal_fontSize <- input$goal_fontSize
      goal_values$goal_mixsize <- input$goal_mixsize
      goal_values$goal_maxsize <- input$goal_maxsize
      
      # browser()
      
      if (input$goal_tools=="Bar chart")
      {
        output$go_plot_scroll_box <- renderUI({
          if (is.null(goal$goal_go)){
            return (NULL)
          }
          else
          {div(style ="overflow-y: scroll; overflow-x: scroll;max-height:700px;border: 1px solid #ddd; padding: 20px;",
               highchartOutput(session$ns("goal_go_plot"),width=goal_values$goal_plotbarWidth,height=goal_values$goal_plotbarHeight)
               #  highchartOutput(session$ns("go_plot"))
               # highchartOutput("go_plot",width="1500px",height="800px")
          )}
        })
        # browser()
        output$goal_go_plot <- renderHighchart({
          my_bar(goal$goal_go,go_term = goal_values$go_term ,color_min = goal_values$goal_col1_bar, color_max = goal_values$goal_col2_bar,
                 fontSize=goal_values$goal_fontSize,sourceWidth=goal_values$goal_plotbarWidth,sourceHeight=goal_values$goal_plotbarHeight)
          # my_bar(go_data())
        })
        
        output$go_download <- renderUI({
          if (is.null(goal$goal_go)) {
            return(NULL)
          }
          div(class = "dropdown",
              tags$button(class = "btn btn-custom dropdown-toggle", type = "button", id = "dropdownMenuButton",
                          `data-bs-toggle` = "dropdown", `aria-haspopup` = "true", `aria-expanded` = "false",
                          "Download  Data ", span(class = "caret")
              ),
              div(class = "dropdown-menu", `aria-labelledby` = "dropdownMenuButton",
                  downloadLink(session$ns("download_goal_go_csv"), "CSV", class = "btn-link"),
                  downloadLink(session$ns("download_goal_go_txt"), "TXT", class = "btn-link")
              ))
          
        })
        
        
        output$go_result <- DT::renderDataTable({
          if (is.null(goal$goal_go)){
            return(NULL)}
          DT::datatable(goal$goal_go,
                        options = list(pageLength =10, scrollX = TRUE))}
        )
        
        output$download_goal_go_csv <- downloadHandler(
          filename = function() {
            paste("result", Sys.Date(), ".csv", sep="")
          },
          content = function(file) {
            # load your_data_frame 
            write.csv(goal$goal_go, file)
          }
        )
        
        # downlaoad txt
        output$download_goal_go_txt <- downloadHandler(
          filename = function() {
            paste("result", Sys.Date(), ".txt", sep="")
          },
          content = function(file) {
            # load your_data_frame
            write.table(goal$goal_go, file, row.names = FALSE)
          }
        )
        
      }
      else if (input$goal_tools=="Bubble char") {
        output$go_plot_scroll_box <- renderUI({
          if (is.null(goal$goal_go)){
            return (NULL)
          }
          else
          {div(style ="overflow-y: scroll; overflow-x: scroll;max-height:700px;border: 1px solid #ddd; padding: 20px;",
               highchartOutput(session$ns("goal_go_plot"),width=goal_values$goal_plotbubbleWidth,height=goal_values$goal_plotbubbleHeight)
               #  highchartOutput(session$ns("go_plot"))
               # highchartOutput("go_plot",width="1500px",height="800px")
          )}
        })
        
        output$goal_go_plot <- renderHighchart({
          my_Bubble_plot(goal$goal_go,go_term = goal_values$go_term,color_min = goal_values$goal_col1_bb, color_max = goal_values$goal_col2_bb,
                         sourceWidth=goal_values$goal_plotbubbleWidth,sourceHeight=goal_values$goal_plotbubbleHeight,
                         fontSize=goal_values$goal_Bubble_fontSize,minsize=goal_values$goal_mixsize,maxsize=goal_values$goal_maxsize)
        })
        
        output$go_download <- renderUI({
          if (is.null(goal$goal_go)) {
            return(NULL)
          }
          div(class = "dropdown",
              tags$button(class = "btn btn-custom dropdown-toggle", type = "button", id = "dropdownMenuButton",
                          `data-bs-toggle` = "dropdown", `aria-haspopup` = "true", `aria-expanded` = "false",
                          "Download  Data ", span(class = "caret")
              ),
              div(class = "dropdown-menu", `aria-labelledby` = "dropdownMenuButton",
                  downloadLink(session$ns("download_goal_go_bubble_csv"), "CSV", class = "btn-link"),
                  downloadLink(session$ns("download_goal_go_bubble_txt"), "TXT", class = "btn-link")
              ))
          
        })
        
        
        output$go_result <- DT::renderDataTable({
          if (is.null(goal$goal_go)){
            return(NULL)}
          DT::datatable(goal$goal_go,
                        options = list(pageLength =10, scrollX = TRUE))}
        )
        
        output$download_goal_go_bubble_csv <- downloadHandler(
          filename = function() {
            paste("result", Sys.Date(), ".csv", sep="")
          },
          content = function(file) {
            # load your_data_frame 
            write.csv(goal$goal_go, file)
          }
        )
        
        # downlaoad txt
        output$download_goal_go_bubble_txt <- downloadHandler(
          filename = function() {
            paste("result", Sys.Date(), ".txt", sep="")
          },
          content = function(file) {
            # load your_data_frame
            write.table(goal$goal_go, file, row.names = FALSE)
          }
        )
        
      }
      
    }
    )
    
    
  }
  #else sentence end part 
  
  )  
  #observer
  #observe
}