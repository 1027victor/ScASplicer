#' dynamics UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
dynamics_UI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      
      fileInput(ns("file"), label = p("File input:",style="color:black; text-align:center")),
      fluidRow(column(4,actionButton(ns("progress"), "calculate",class='btn-info')),
               column(4,actionButton(ns("action"), "Visualization",class='btn- primary',icon=icon('square-poll-horizontal'))))
      
      
    ),
    
    mainPanel(splitLayout(div(helpText("Show the Modality dynamics results in this view.",style="font-weight:bold;font-size:24px;color:#303133 !important;font-family: OpenSans-Semibold, Helvetica Neue, Helvetica, Arial, sans-serif !important;"), 
                              uiOutput(ns("modality_scroll_box")),
                              br(),
                              uiOutput(ns("modality_download")),
                              br(),
                              uiOutput(ns("modality_assign_dynamics")),
                              # fluidRow(
                              # column(6,DT::dataTableOutput(ns('modality_table'))),
                              # column(6,plotOutput(ns("modality_plot")))),
                              helpText("Show the Gene-splicing dynamics in this view.",style="font-weight:bold;font-size:24px;color:#303133 !important;font-family: OpenSans-Semibold, Helvetica Neue, Helvetica, Arial, sans-serif !important;"),
                              uiOutput(ns("Gene_splicing_scroll_box")),
                              br(),
                              uiOutput(ns("Gene_splicing_download")),
                              br(),
                              uiOutput(ns('Gene_splicing_dynamics')),
                              
    ),  wellPanel(
      div(colourpicker::colourInput(ns("col1"),"Select color","#1B9E77"),
          colourpicker::colourInput(ns("col2"),"Select color","#D95F02"),
          colourpicker::colourInput(ns("col3"),"Select color","#7570B3"),
          colourpicker::colourInput(ns("col4"),"Select color","#E7298A"),
          # colourpicker::colourInput("col7","Select color","#A6761D"),
          wellPanel(div(p("Pie Plot(Assign dynamics)",style="font-weight: bold; color: #000000; font-size: 16px;"),
                        fluidRow(column(6,textInput(ns("plot1width"),"Plot width",value = 400)),
                                 column(6,textInput(ns("plot1height"),"Plot height",value = 400))),
                        p('Volcano Plot(Assign dynamics)',style="font-weight: bold; color: #000000; font-size: 16px;"),
                        fluidRow(column(6,textInput(ns("plotwidth"),"Plot width",value = 6)),
                                 column(6,textInput(ns("plotheight"),"Plot height",value = 5))),
                        radioButtons(ns('extPlot'), p('Plot output format',style="font-weight: bold; color: #000000; font-size: 16px;"),choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                        downloadButton(ns("plot1downloadData"),"Download Plot"))),
          br(),
          wellPanel(div(p("Pie Plot(Gene_splicing dynamics)",style="font-weight: bold; color: #000000; font-size: 16px;"),
                        fluidRow(column(6,textInput(ns("plot1gpwidth"),"Plot width",value = 400)),
                                 column(6,textInput(ns("plot1gpheight"),"Plot height",value = 400))),
                        p('Volcano Plot(Gene_splicing dynamics)',style="font-weight: bold; color: #000000; font-size: 16px;"),
                        fluidRow(column(6,textInput(ns("plotgpwidth"),"Plot width",value = 12)),
                                 column(6,textInput(ns("plotgpheight"),"Plot height",value = 8))),
                        radioButtons(ns('extgpPlot'), p('Plot output format',style="font-weight: bold; color: #000000; font-size: 16px;"),choices = c("PNG"='png', 'PDF'='pdf','JPEG'='jpeg'), inline = T),
                        downloadButton(ns("plot1gpdownloadData"),"Download Plot"))))),
    
    cellWidths = c("70%","30%"))
    )
  )
}

dynamics_mainServer <- function(id) {
  moduleServer(id, function(input, output,session) {
    
    data_modality_dynamics <- reactiveValues()
    
    # process data
    observeEvent(input$progress,{
      withProgress(message = 'processing...', value = 0, {
        infile <- input$file
        if (is.null(infile)) {
          return(NULL)
        }
        load(infile$datapath)
        data<-list()
        data_isoswitch<-list()
        marvel_object<-list()
        marvel_object_isoswitch<-list()
        merge_data<-NULL
        merge_data_isoswitch<-NULL
        l<-1
        df.pheno <- marvel$SplicePheno
        cell_meta_type<-as.data.frame(table(df.pheno$cell.type))
        cell_types<-unique(df.pheno$cell.type)
        total_iterations <- (length(cell_types) - 1) * length(cell_types) / 2
        progress <- 0
        # Differential analysis
        for (i in 1:(length(cell_types) - 1)) {
          num_one<-(cell_meta_type[cell_meta_type$Var1==cell_types[i],]$Freq)
          for (j in (i + 1):length(cell_types)) {
            num_two<-(cell_meta_type[cell_meta_type$Var1==cell_types[j],]$Freq)
            cell.group.g1 <- df.pheno[which(df.pheno$cell.type==cell_types[i]), "sample.id"]
            cell.group.g2 <- df.pheno[which(df.pheno$cell.type==cell_types[j]), "sample.id"]
            marvel <- CompareValues(MarvelObject=marvel,
                                    cell.group.g1=cell.group.g1,
                                    cell.group.g2=cell.group.g2,
                                    min.cells=3,
                                    method="wilcox",
                                    method.adjust="fdr",
                                    level="gene",
                                    show.progress=FALSE
            )
            
            
            # Differential splicing analysis
            num<-round(min(num_one,num_two)/4)+1
            if(num>25){num=25}
            marvel <- CompareValues(MarvelObject=marvel,
                                    cell.group.g1=cell.group.g1,
                                    cell.group.g2=cell.group.g2,
                                    min.cells=num,
                                    method=c("ad", "dts"),
                                    method.adjust="fdr",
                                    level="splicing",
                                    event.type=c("SE", "MXE", "RI", "A5SS", "A3SS"),
                                    show.progress=FALSE
            )
            #Differential (spliced) gene analysis
            marvel <- CompareValues(MarvelObject=marvel,
                                    cell.group.g1=cell.group.g1,
                                    cell.group.g2=cell.group.g2,
                                    psi.method=c("ad", "dts"),
                                    #set shiny parameter
                                    psi.pval=c(0.10, 0.10),
                                    #set shiny parameter
                                    psi.delta=0,
                                    method.de.gene="wilcox",
                                    method.adjust.de.gene="fdr",
                                    downsample=FALSE,
                                    show.progress=FALSE,
                                    level="gene.spliced"
            )
            # browser()
            #Modality_dynamics
            marvel <- ModalityChange(MarvelObject=marvel,
                                     method=c("ad", "dts"),
                                     psi.pval=c(0.10, 0.10)
            )
            
            
            data[[l]]<-marvel$DE$Modality$Plot.Stats
            #name
            marvel_object[[l]]<-marvel
            data_modal<-marvel$DE$Modality$Table
            data_modal$group_one<-paste(cell_types[i])
            data_modal$group_two<-paste(cell_types[j])
            names(data)[l]<-paste(cell_types[i], " vs ", cell_types[j])
            #name
            names(marvel_object)[l]<-paste(cell_types[i], " vs ", cell_types[j])
            
            
            
            #Assign dynamics
            # browser()
            marvel <- IsoSwitch(MarvelObject=marvel,
                                method=c("ad", "dts"),
                                psi.pval=c(0.10, 0.10),
                                psi.delta=0,
                                gene.pval=0.10,
                                gene.log2fc=0.5
            )
            
            data_isoswitch[[l]]<-marvel$DE$Cor$Plot.Stats
            marvel_object_isoswitch[[l]]<-marvel
            data_gene_splicing<-marvel$DE$Cor$Table
            data_gene_splicing$group_one<-paste(cell_types[i])
            data_gene_splicing$group_two<-paste(cell_types[j])
            names(data_isoswitch)[l]<-paste(cell_types[i], " vs ", cell_types[j])
            #name
            names(marvel_object_isoswitch)[l]<-paste(cell_types[i], " vs ", cell_types[j])
            
            
            l<-l+1
            if(is.null(merge_data)) {
              merge_data <- data_modal
            } else {
              merge_data<- rbind(merge_data,data_modal)
            }
            
            if(is.null(merge_data_isoswitch)) {
              merge_data_isoswitch <- data_gene_splicing
            } else {
              merge_data_isoswitch<- rbind(merge_data_isoswitch,data_gene_splicing)
            }
            
            progress <- progress + 1
            setProgress(value = progress / total_iterations,message = sprintf("Processing progress: %d%%", round((progress / total_iterations) * 100)))
          }
        }
        # showNotification("complete! ",closeButton = TRUE,
        #     duration = 60,type = c("message"))
        showModal(modalDialog(
          title = "",
          "The data processing is complete!",
          footer = modalButton("close")
        ))
        
      })
      data_modality_dynamics$modal=data
      data_modality_dynamics$modal_merge=merge_data
      data_modality_dynamics$MARVEL=marvel_object
      
      data_modality_dynamics$isoswitch=data_isoswitch
      data_modality_dynamics$isoswitch_merge=merge_data_isoswitch
      data_modality_dynamics$MARVEL_isoswitch=marvel_object_isoswitch
      
    },ignoreInit = TRUE)
    
    
    
    
    radarplot_modality<-reactive({
      req(length(data_modality_dynamics$modal) > 0)
      chart_outputs <- lapply(seq_along(data_modality_dynamics$modal), function(i) {
        chart_id <- paste("chart_modality", i, sep = "_")
        withSpinner(highchartOutput(session$ns(chart_id), height = '400px'),type=8)
      })
      
      chart_rows <- lapply(seq(1, length(chart_outputs), by = 2), function(i) {
        fluidRow(
          column(6,chart_outputs[[i]]),
          if (i + 1 <= length(chart_outputs)) column(6, chart_outputs[[i + 1]])
        )
      })
      # function
      do.call(tagList, chart_rows)
    })
    
    
    #IsoSwitch
    
    radarplot_IsoSwitch<-reactive({
      req(length(data_modality_dynamics$isoswitch) > 0)
      chart_outputs <- lapply(seq_along(data_modality_dynamics$isoswitch), function(i) {
        chart_id <- paste("chart_isoswitch", i, sep = "_")
        withSpinner(highchartOutput(session$ns(chart_id), height = '400px'),type=8)
      })
      
      chart_rows <- lapply(seq(1, length(chart_outputs), by = 2), function(i) {
        fluidRow(
          column(6,chart_outputs[[i]]),
          if (i + 1 <= length(chart_outputs)) column(6, chart_outputs[[i + 1]])
        )
      })
      do.call(tagList, chart_rows)
    })
    
    
    
    
    observeEvent(input$action,{
      warm_colors <- c(input$col1,input$col2,input$col3,input$col4)
      lapply(seq_along(data_modality_dynamics$modal), function(i) {
        chart_id <- paste("chart_modality", i, sep = "_")
        data <- data_modality_dynamics$modal[[i]]
        output[[chart_id]] <- renderHighchart({
          
          
          highchart() %>%
            hc_chart(backgroundColor = "#FFFFFF",type = "pie") %>%
            hc_plotOptions(pie = list(innerSize = '60%')) %>%
            hc_colors(warm_colors) %>% 
            hc_series(list(
              name = 'freq',
              data = data %>% mutate(name = modality.change, y =freq) %>% select(name, y) %>% list_parse()
            )) %>%
            hc_exporting(enabled = TRUE,sourceWidth = input$plot1width, sourceHeight = input$plot1height) %>%
            
            # hc_size(height = 600, width = 800)%>%
            hc_title(text = names(data_modality_dynamics$modal)[i])
        })
      })
      
      output$modality_scroll_box <- renderUI({
        if (length(data_modality_dynamics$modal) > 0) {
          div(style ="overflow-y: scroll; max-height: 800px; border: 1px solid #ddd; padding: 20px;",
              uiOutput(session$ns("charts"))
          )
        }
      })
      output$charts<-renderUI({radarplot_modality()})
      
      
      # download data
      
      output$modality_download <- renderUI({
        if (is.null(data_modality_dynamics$modal_merge)) {
          return(NULL)
        }
        div(class = "dropdown",
            tags$button(class = "btn btn-custom dropdown-toggle", type = "button", id = "dropdownMenuButton",
                        `data-bs-toggle` = "dropdown", `aria-haspopup` = "true", `aria-expanded` = "false",
                        "Download  Data ", span(class = "caret")
            ),
            div(class = "dropdown-menu", `aria-labelledby` = "dropdownMenuButton",
                downloadLink(session$ns("download_modality_csv"), "CSV", class = "btn-link"),
                downloadLink(session$ns("download_modality_txt"), "TXT", class = "btn-link")
            ))
        
      })
      
      
      output$modality_assign_dynamics<-renderUI({
        div(DT::dataTableOutput(session$ns('modality_table')),
            br(),
            plotOutput(session$ns("modality_plot"),width='800px',height = "300px"))
      })
      
      #select data plot
      
      
      output$modality_table <- DT::renderDataTable({
        if (is.null(data_modality_dynamics$modal_merge)){
          return(NULL)}
        DT::datatable(data_modality_dynamics$modal_merge,selection = 'single',
                      options = list(pageLength =10, scrollX = TRUE))}
      )
      
      
      # start plot
      plot<-reactive({
        selectedRow <- input$modality_table_rows_selected
        if (length(selectedRow)) {
          # tran_id
          selectedData <- data_modality_dynamics$modal_merge$tran_id[selectedRow]
          name_one <- data_modality_dynamics$modal_merge$group_one[selectedRow]
          name_two <- data_modality_dynamics$modal_merge$group_two[selectedRow]
          # cell.group.list <- setNames(list(sample.ids.1, sample.ids.2), c(name_one, name_two))
          index<-paste(name_one, " vs ",  name_two)
          df.pheno <- (data_modality_dynamics$MARVEL)[[index]]$SplicePheno
          cell_meta_type<-as.data.frame(table(df.pheno$cell.type))
          num_one<-(cell_meta_type[cell_meta_type$Var1==name_one,]$Freq)
          num_two<-(cell_meta_type[cell_meta_type$Var1==name_two,]$Freq)
          sample.ids.1 <- df.pheno[which(df.pheno$cell.type==name_one), "sample.id"]
          sample.ids.2 <- df.pheno[which(df.pheno$cell.type==name_two), "sample.id"]
          cell.group.list <- setNames(list(sample.ids.1, sample.ids.2), c(name_one, name_two))
          min_cell_num <- round(min(num_one,num_two)/4)+1
          if(min_cell_num>25)
          {min_cell_num=25}
          # name_one <- data_modality_dynamics$modal_merge$group_one[selectedRow]
          # name_two <- data_modality_dynamics$modal_merge$group_two[selectedRow]
          # cell.group.list <- setNames(list(sample.ids.1, sample.ids.2), c(name_one, name_two))
          # index<-paste(name_one, " vs ",  name_two)
          marvel <- PlotValues(MarvelObject=(data_modality_dynamics$MARVEL)[[index]],
                               cell.group.list=cell.group.list,
                               feature=selectedData,
                               xlabels.size=12,
                               level="splicing",
                               min.cells=min_cell_num,
          )
          marvel$adhocPlot$PSI
        }
        else {
          return(NULL)
        }
      })
      
      output$modality_plot <- renderPlot({
        if (is.null(data_modality_dynamics$modal_merge)) {
          return(NULL)
        }
        # selectedRow <- input$modality_table_rows_selected
        # if (length(selectedRow)) {
        #   selectedData <- data_modality_dynamics$modal_merge$tran_id[selectedRow]
        #   df.pheno <- data_modality_dynamics$MARVEL$SplicePheno
        #   cell_meta_type<-as.data.frame(table(df.pheno$cell.type))
        #   num_one<-(cell_meta_type[cell_meta_type$Var1==data_modality_dynamics$modal_merge$group_one[selectedRow],]$Freq)
        #   num_two<-(cell_meta_type[cell_meta_type$Var1==data_modality_dynamics$modal_merge$group_two[selectedRow],]$Freq)
        #   sample.ids.1 <- df.pheno[which(df.pheno$cell.type==data_modality_dynamics$modal_merge$group_one[selectedRow]), "sample.id"]
        #   sample.ids.2 <- df.pheno[which(df.pheno$cell.type==data_modality_dynamics$modal_merge$group_two[selectedRow]), "sample.id"]
        #   min_cell_num <- round(min(num_one,num_two)/4)+1
        #   if(min_cell_num>25)
        #   {min_cell_num=25}
        #   name_one <- data_modality_dynamics$modal_merge$group_one[selectedRow]
        #   name_two <- data_modality_dynamics$modal_merge$group_two[selectedRow]
        #   cell.group.list <- setNames(list(sample.ids.1, sample.ids.2), c(name_one, name_two))
        # 
        #   marvel <- PlotValues(MarvelObject=data_modality_dynamics$MARVEL,
        #                cell.group.list=cell.group.list,
        #                feature=selectedData,
        #                xlabels.size=12,
        #                level="splicing",
        #                min.cells=min_cell_num,
        #                )
        #   marvel$adhocPlot$PSI
        plot()
      })
      
      
      output$download_modality_csv <- downloadHandler(
        filename = function() {
          paste("result", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          # load your_data_frame 
          write.csv(data_modality_dynamics$modal_merge, file)
        }
      )
      
      # download txt
      output$download_modality_txt <- downloadHandler(
        filename = function() {
          paste("result", Sys.Date(), ".txt", sep="")
        },
        content = function(file) {
          # load your_data_frame
          write.table(data_modality_dynamics$modal_merge, file, row.names = FALSE)
        }
      )
      
      
      # download plot
      output$plot1downloadData<-downloadHandler(
        filename = function() {
          paste('Modality_dynamics', Sys.Date(), '.',input$extPlot, sep='')
        },
        content=function(file){
          ggsave(filename = file,plot = plot(), width = as.numeric(input$plotwidth), height = as.numeric(input$plotheight), dpi = 300)
        })
      
      
      
      
      # Gene_splicing dynamics
      
      
      lapply(seq_along(data_modality_dynamics$isoswitch), function(i) {
        chart_id <- paste("chart_isoswitch", i, sep = "_")
        data <- data_modality_dynamics$isoswitch[[i]]
        output[[chart_id]] <- renderHighchart({
          
          
          highchart() %>%
            hc_chart(backgroundColor = "#FFFFFF",type = "pie") %>%
            hc_plotOptions(pie = list(innerSize = '60%')) %>%
            hc_colors(warm_colors) %>% 
            hc_series(list(
              name = 'freq',
              data = data %>% mutate(name = cor, y =freq) %>% select(name, y) %>% list_parse()
            )) %>%
            hc_exporting(enabled = TRUE,sourceWidth = input$plot1gpwidth, sourceHeight = input$plot1gpheight) %>%
            
            # hc_size(height = 600, width = 800)%>%
            hc_title(text = names(data_modality_dynamics$isoswitch)[i])
        })
      })
      
      output$Gene_splicing_scroll_box <- renderUI({
        if (length(data_modality_dynamics$isoswitch) > 0) {
          div(style ="overflow-y: scroll; max-height: 800px; border: 1px solid #ddd; padding: 20px;",
              uiOutput(session$ns("charts_isoswitch"))
          )
        }
      })
      output$charts_isoswitch<-renderUI({radarplot_IsoSwitch()})
      
      
      # download data
      
      output$Gene_splicing_download <- renderUI({
        if (is.null(data_modality_dynamics$isoswitch_merge)) {
          return(NULL)
        }
        div(class = "dropdown",
            tags$button(class = "btn btn-custom dropdown-toggle", type = "button", id = "dropdownMenuButton",
                        `data-bs-toggle` = "dropdown", `aria-haspopup` = "true", `aria-expanded` = "false",
                        "Download  Data ", span(class = "caret")
            ),
            div(class = "dropdown-menu", `aria-labelledby` = "dropdownMenuButton",
                downloadLink(session$ns("download_isoswitch_csv"), "CSV", class = "btn-link"),
                downloadLink(session$ns("download_isoswitch_txt"), "TXT", class = "btn-link")
            ))
        
      })
      
      
      output$Gene_splicing_dynamics<-renderUI({
        div(DT::dataTableOutput(session$ns('Gene_table')),
            br(),
            DT::dataTableOutput(session$ns('splicing_table')),
            br(),
            plotOutput(session$ns("Gene_splicing_plot"),height = "600px"))
      })
      
      #select data plot
      
      
      output$Gene_table <- DT::renderDataTable({
        if (is.null(data_modality_dynamics$isoswitch_merge)){
          return(NULL)}
        DT::datatable(data_modality_dynamics$isoswitch_merge,selection = 'single',
                      options = list(pageLength =10, scrollX = TRUE))}
      )
      
      
      output$splicing_table <- DT::renderDataTable({
        if (is.null(data_modality_dynamics$modal_merge)){
          return(NULL)}
        DT::datatable(data_modality_dynamics$modal_merge,selection = 'single',
                      options = list(pageLength =10, scrollX = TRUE))}
      )
      
      
      # start plot
      isoswitch_plot<-reactive({
        #同时左右选中两行信息
        selectedRow <- input$splicing_table_rows_selected
        selectedRow_gene<-input$Gene_table_rows_selected
        #修改if条件，同时满足两个条件
        if (length(selectedRow)>0 && length(selectedRow_gene)>0 ) {
          # browser()
          selectedData_gene <- data_modality_dynamics$isoswitch_merge$gene_short_name[selectedRow_gene]
          selectedData <- data_modality_dynamics$modal_merge$tran_id[selectedRow]
          name_one <- data_modality_dynamics$isoswitch_merge$group_one[selectedRow]
          name_two <- data_modality_dynamics$isoswitch_merge$group_two[selectedRow]
          index<-paste(name_one, " vs ",  name_two)
          df.pheno <- (data_modality_dynamics$MARVEL_isoswitch)[[index]]$SplicePheno
          cell_meta_type<-as.data.frame(table(df.pheno$cell.type))
          num_one<-(cell_meta_type[cell_meta_type$Var1==name_one,]$Freq)
          num_two<-(cell_meta_type[cell_meta_type$Var1==name_two,]$Freq)
          sample.ids.1 <- df.pheno[which(df.pheno$cell.type==name_one), "sample.id"]
          sample.ids.2 <- df.pheno[which(df.pheno$cell.type==name_two), "sample.id"]
          cell.group.list <- setNames(list(sample.ids.1, sample.ids.2), c(name_one, name_two))
          
          min_cell_num <- round(min(num_one,num_two)/4)+1
          if(min_cell_num>25)
          {min_cell_num=25}
          # name_one <- data_modality_dynamics$modal_merge$group_one[selectedRow]
          # name_two <- data_modality_dynamics$modal_merge$group_two[selectedRow]
          # cell.group.list <- setNames(list(sample.ids.1, sample.ids.2), c(name_one, name_two))
          # index<-paste(name_one, " vs ",  name_two)
          df.feature <- (data_modality_dynamics$MARVEL_isoswitch)[[index]]$GeneFeature
          gene_id <- df.feature[which(df.feature$gene_short_name==selectedData_gene), "gene_id"]
          marvel <- PlotValues(MarvelObject=(data_modality_dynamics$MARVEL_isoswitch)[[index]],
                               cell.group.list=cell.group.list,
                               feature=gene_id,
                               maintitle="gene_short_name",
                               xlabels.size=12,
                               level="gene"
          )
          
          plot.1_gene <- marvel$adhocPlot$Exp
          
          marvel <- PlotValues(MarvelObject=(data_modality_dynamics$MARVEL)[[index]],
                               cell.group.list=cell.group.list,
                               feature=selectedData,
                               xlabels.size=12,
                               level="splicing",
                               min.cells=min_cell_num,
          )
          
          plot.1_splicing <- marvel$adhocPlot$PSI
          grid.arrange(plot.1_gene, plot.1_splicing,
                       nrow=1)
        }
        else {
          return(NULL)
        }
      })
      
      output$Gene_splicing_plot <- renderPlot({
        if (is.null(data_modality_dynamics$isoswitch_merge)) {
          return(NULL)
        }
        # selectedRow <- input$modality_table_rows_selected
        # if (length(selectedRow)) {
        #   selectedData <- data_modality_dynamics$modal_merge$tran_id[selectedRow]
        #   df.pheno <- data_modality_dynamics$MARVEL$SplicePheno
        #   cell_meta_type<-as.data.frame(table(df.pheno$cell.type))
        #   num_one<-(cell_meta_type[cell_meta_type$Var1==data_modality_dynamics$modal_merge$group_one[selectedRow],]$Freq)
        #   num_two<-(cell_meta_type[cell_meta_type$Var1==data_modality_dynamics$modal_merge$group_two[selectedRow],]$Freq)
        #   sample.ids.1 <- df.pheno[which(df.pheno$cell.type==data_modality_dynamics$modal_merge$group_one[selectedRow]), "sample.id"]
        #   sample.ids.2 <- df.pheno[which(df.pheno$cell.type==data_modality_dynamics$modal_merge$group_two[selectedRow]), "sample.id"]
        #   min_cell_num <- round(min(num_one,num_two)/4)+1
        #   if(min_cell_num>25)
        #   {min_cell_num=25}
        #   name_one <- data_modality_dynamics$modal_merge$group_one[selectedRow]
        #   name_two <- data_modality_dynamics$modal_merge$group_two[selectedRow]
        #   cell.group.list <- setNames(list(sample.ids.1, sample.ids.2), c(name_one, name_two))
        # 
        #   marvel <- PlotValues(MarvelObject=data_modality_dynamics$MARVEL,
        #                cell.group.list=cell.group.list,
        #                feature=selectedData,
        #                xlabels.size=12,
        #                level="splicing",
        #                min.cells=min_cell_num,
        #                )
        #   marvel$adhocPlot$PSI
        isoswitch_plot()
      })
      
      
      output$download_isoswitch_csv <- downloadHandler(
        filename = function() {
          paste("result", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          # load your_data_frame 
          write.csv(data_modality_dynamics$isoswitch_merge, file)
        }
      )
      
      # download txt
      output$download_isoswitch_txt <- downloadHandler(
        filename = function() {
          paste("result", Sys.Date(), ".txt", sep="")
        },
        content = function(file) {
          # load your_data_frame
          write.table(data_modality_dynamics$isoswitch_merge, file, row.names = FALSE)
        }
      )
      
      
      # download plot
      output$plot1gpdownloadData<-downloadHandler(
        filename = function() {
          paste('Gene_splicing_dynamics', Sys.Date(), '.',input$extgpPlot, sep='')
        },
        content=function(file){
          ggsave(filename = file,plot = isoswitch_plot(), width = as.numeric(input$plotgpwidth), height = as.numeric(input$plotgpheight), dpi = 300)
        })
      
      
    })
    
    
    
    
    
  }
  )  
}