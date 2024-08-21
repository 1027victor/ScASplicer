#' modality_analysis UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
modal_UI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      
      
      
      fileInput(ns("modal_file"), label = p("File input:",style="color:black; text-align:center")),
      fluidRow( column(4,actionButton(ns("modal_progress"), "calculate",class='btn-info')),
                column(4,actionButton(ns("modal_action"), "Visualization",class='btn-primary',icon=icon('square-poll-horizontal'))))
      
      
    ),
    mainPanel(splitLayout(div(helpText("Show the results in this view.",style="font-weight:bold;font-size:24px;color:#303133 !important;font-family: OpenSans-Semibold, Helvetica Neue, Helvetica, Arial, sans-serif !important;"),
                              uiOutput(ns("modal_conditional_scroll_box")),
                              # c(rep(list(br()), 10)),
                              
                              
                              helpText("Show the results in this view.",style="font-weight:bold;font-size:24px;color:#303133 !important;font-family: OpenSans-Semibold, Helvetica Neue, Helvetica, Arial, sans-serif !important;"), uiOutput(ns("modal_bar_conditional_scroll_box")),
                              helpText("Show the output_data in this view.",style="font-weight:bold;font-size:24px;color:#303133 !important;font-family: OpenSans-Semibold, Helvetica Neue, Helvetica, Arial, sans-serif !important;"),
                              br(),
                              uiOutput(ns("modal")),
                              br(),
                              DT::dataTableOutput(ns('modal_result'))),  wellPanel(div(colourpicker::colourInput(ns("colm1"),"Select color","#1B9E77"),
                                                                                                                                            colourpicker::colourInput(ns("colm2"),"Select color","#D95F02"),
                                                                                                                                            colourpicker::colourInput(ns("colm3"),"Select color","#7570B3"),
                                                                                                                                            colourpicker::colourInput(ns("colm4"),"Select color","#E7298A"),
                                                                                                                                            colourpicker::colourInput(ns("colm5"),"Select color","#66A61E"),
                                                                                                                                            colourpicker::colourInput(ns("colm6"),"Select color","#E6AB02"),
                                                                                                                                            colourpicker::colourInput(ns("colm7"),"Select color","#A6761D"),
                                                                                                                                            fluidRow(column(6,textInput(ns("plot1width_modal"),"Plot width",value = 400)),
                                                                                                                                                     column(6,textInput(ns("plot1height_modal"),"Plot height",value = 400))))),
                          
                          cellWidths = c("70%","30%"))
    )
  )
}


modal_mainServer <- function(id) {
  moduleServer(id, function(input, output,session) {
    data_list_modal<-reactiveValues()
    observeEvent(input$modal_progress,{
      withProgress(message = 'processing...', value = 0, {
        infile <- input$modal_file
        if (is.null(infile)) {
          return(NULL)
        }
        load(infile$datapath)
        df.pheno <- marvel$SplicePheno
        cell_meta_type<-as.data.frame(table(df.pheno$cell.type))
        # new add  
        merge_data<-NULL
        #
        
        data_modal_pie<-list()
        data_modal_bar<-list()
        j<-1
        # Iterate over unique cell types
        total_iterations <- length(unique(df.pheno$cell.type))
        progress <- 0
        # withProgress(message = 'processing...', value = 0, {
        for (i in sort(unique(df.pheno$cell.type))){
          
          # Define sample ids
          sample.ids <- df.pheno[df.pheno$cell.type == i, "sample.id"]
          num<-round((cell_meta_type[cell_meta_type$Var1==i,]$Freq)/4)+1
          if (num>25){num=25}
          marvel <- AssignModality(MarvelObject=marvel,
                                   sample.ids=sample.ids,
                                   min.cells=num,
                                   seed=1
          )
          
          marvel <- PropModality(MarvelObject=marvel,
                                 modality.column="modality.bimodal.adj",
                                 modality.type="extended",
                                 event.type=c("SE", "MXE", "RI", "A5SS", "A3SS"),
                                 across.event.type=FALSE
          )
          marvel <- PropModality(MarvelObject=marvel,
                                 modality.column="modality.bimodal.adj",
                                 modality.type="extended",
                                 event.type=c("SE", "MXE", "RI", "A5SS", "A3SS"),
                                 across.event.type=TRUE,
                                 prop.test="chisq",
                                 prop.adj="fdr",
                                 xlabels.size=8
          )
          # new add
          data_modality_result<-marvel$Modality$Results[,c("tran_id", "event_type", "gene_id", "gene_short_name", "modality.bimodal.adj")]
          data_modality_result$cell_type<-i
          #
          
          data_modal_pie[[j]]<-marvel$Modality$Prop$DoughnutChart$Table
          names(data_modal_pie)[j]<-i
          data_modal_bar[[j]] <- marvel$Modality$Prop$BarChart$Table
          names(data_modal_bar)[j]<-i
          j<-j+1
          
          
          # new add 
          
          if(is.null(merge_data)) {
            merge_data <- data_modality_result
          } else {
            merge_data <- rbind(merge_data,data_modality_result)
          }
          #
          
          progress <- progress + 1
          setProgress(value = progress / total_iterations,message = sprintf("Processing progress: %d%%", round((progress / total_iterations) * 100)))
        }
        # showNotification("complete! ", type = "warning")
        showModal(modalDialog(
          title = "",
          "The data processing is complete!",
          footer = modalButton("close")
        ))
      }
      )
      # list(data_one=data_modal_pie,data_two=data_modal_bar)
      data_list_modal$data_one=data_modal_pie
      data_list_modal$data_two=data_modal_bar
      # new add
      data_list_modal$result=merge_data
      #
    },ignoreInit = TRUE)
    
    radarplot_modal<-reactive({
      req(length(data_list_modal$data_one) > 0)
      chart_outputs <- lapply(seq_along(data_list_modal$data_one), function(i) {
        chart_id <- paste("chart_modal", i, sep = "_")
        highchartOutput(session$ns(chart_id), height = "400px")
      })
      
      chart_rows <- lapply(seq(1, length(chart_outputs), by = 2), function(i) {
        fluidRow(
          column(6,chart_outputs[[i]]),
          if (i + 1 <= length(chart_outputs)) column(6, chart_outputs[[i + 1]])
          # if (i + 1 <= length(chart_outputs)) column(4, chart_outputs[[i + 2]])
        )
      })
      do.call(tagList, chart_rows)
    })
    
    
    observeEvent(input$modal_action,{
      warm_colors <- c(input$colm1,input$colm2,input$colm3,input$colm4,input$colm5,input$colm6,input$colm7)
      lapply(seq_along(data_list_modal$data_one), function(i) {
        chart_id <- paste("chart_modal", i, sep = "_")
        data <- data_list_modal$data_one[[i]]
        output[[chart_id]] <- renderHighchart({
          
          
          highchart() %>%
            hc_chart(backgroundColor = "#FFFFFF",type = "pie") %>%
            hc_plotOptions(pie = list(innerSize = '60%')) %>%
            hc_colors(warm_colors) %>%
            hc_series(list(
              name = 'freq',
              data = data %>% mutate(name = data$modality, y =data$freq) %>% select(name, y) %>% list_parse()
            )) %>%
            hc_exporting(enabled = TRUE) %>%
            hc_title(text = names(data_list_modal$data_one)[i])
        })
      })
      output$modal_conditional_scroll_box <- renderUI({
        if (length(data_list_modal) > 0) {
          div(style ="overflow-y: scroll; max-height: 800px; border: 1px solid #ddd; padding: 20px;",
              uiOutput(session$ns("charts_modal"))
          )
        }
      })
      output$charts_modal<-renderUI({radarplot_modal()})
      
    })
    
    
    
    
    
    
    # modal bar plot
    
    radarplot_modal_bar<-reactive({
      req(length(data_list_modal$data_two) > 0)
      chart_outputs <- lapply(seq_along(data_list_modal$data_two), function(i) {
        chart_id <- paste("chart_modal_bar", i, sep = "_")
        highchartOutput(session$ns(chart_id), height = "400px")
      })
      
      chart_rows <- lapply(seq(1, length(chart_outputs), by = 2), function(i) {
        fluidRow(
          column(6,chart_outputs[[i]]),
          if (i + 1 <= length(chart_outputs)) column(6, chart_outputs[[i + 1]])
          # if (i + 1 <= length(chart_outputs)) column(4, chart_outputs[[i + 2]])
        )
      })
      do.call(tagList, chart_rows)
    })
    
    
    observeEvent(input$modal_action,{
      warm_colors <- c(input$colm1,input$colm2,input$colm3,input$colm4,input$colm5,input$colm6,input$colm7)
      lapply(seq_along(data_list_modal$data_two), function(i) {
        chart_id <- paste("chart_modal_bar", i, sep = "_")
        data <- data_list_modal$data_two[[i]]
        output[[chart_id]] <- renderHighchart({
          
          
          highchart() %>%
            hc_chart(type = "column", backgroundColor = "#FFFFFF") %>%
            hc_title(text = names(data_list_modal$data_two)[i]) %>%
            hc_xAxis(
              categories = unique(data$modality),
              title = list(text = ""),
              lineColor = "black",
              lineWidth = 1,
              gridLineWidth = 0,
              labels = list(style = list(fontSize = "14px", color = "black", fontWeight = "bold"))
            ) %>%
            hc_yAxis(
              title = list(
                text = "%",
                style = list(fontSize = "14px", color = "black", fontWeight = "bold")
              ),
              lineColor = "black",
              lineWidth = 1,
              gridLineWidth = 0,
              labels = list(
                style = list(fontSize = "14px", color = "black", fontWeight = "bold")
              )
            ) %>%
            hc_plotOptions(column = list(pointPadding = 0,borderWidth = 0)) %>%
            hc_add_series(
              data = data,
              type = "column",
              hcaes(x = modality, y = pct, group = event_type),
              showInLegend = TRUE
            ) %>%
            hc_colors(warm_colors) %>%
            hc_legend(title = list(text = "Event Type")) %>%
            hc_exporting(enabled = TRUE,sourceWidth = input$plot1width_modal, sourceHeight = input$plot1height_modal) %>%
            hc_tooltip(shared = FALSE, headerFormat = "", pointFormat = "<b>{point.y:.2f}%</b>")
        })
      })
      
      output$modal_bar_conditional_scroll_box <- renderUI({
        if (length(data_list_modal) > 0) {
          div(style ="overflow-y: scroll; max-height: 800px; border: 1px solid #ddd; padding: 20px;",
              uiOutput(session$ns("charts_modal_bar"))
          )
        }
      })
      output$charts_modal_bar<-renderUI({radarplot_modal_bar()})
      
      output$modal <- renderUI({
        if (is.null(data_list_modal)) {
          return(NULL)
        }
        div(class = "dropdown",
            tags$button(class = "btn btn-custom dropdown-toggle", type = "button", id = "dropdownMenuButton",
                        `data-bs-toggle` = "dropdown", `aria-haspopup` = "true", `aria-expanded` = "false",
                        "Download  Data ", span(class = "caret")
            ),
            div(class = "dropdown-menu", `aria-labelledby` = "dropdownMenuButton",
                downloadLink(session$ns("download_modal_csv"), "CSV", class = "btn-link"),
                downloadLink(session$ns("download_modal_txt"), "TXT", class = "btn-link")
            ))
        
      })
      
      
      output$modal_result <- DT::renderDataTable({
        if (is.null(data_list_modal)){
          return(NULL)}
        DT::datatable(data_list_modal$result,
                      options = list(pageLength =10, scrollX = TRUE))}
      )
      
      output$download_modal_csv <- downloadHandler(
        filename = function() {
          paste("result", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          # load your_data_frame 
          write.csv(data_list_modal$result, file)
        }
      )
      
      # downlaoad txt
      output$download_modal_txt <- downloadHandler(
        filename = function() {
          paste("result", Sys.Date(), ".txt", sep="")
        },
        content = function(file) {
          # load your_data_frame
          write.table(data_list_modal$result, file, row.names = FALSE)
        }
      )
      
      
    
      
      })
    
  }
  )
}