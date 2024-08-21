#' overview UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
overview_UI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      
      
      
      fileInput(ns("file"), label = p("File input:",style="color:black; text-align:center")),
      fluidRow( column(4,actionButton(ns("over_progress"), "calculate",class='btn-info')),
                column(4,actionButton(ns("over_action"), "Visualization",class='btn-primary',icon=icon('square-poll-horizontal'))))
    ),
    
    mainPanel(splitLayout(div(helpText("Show the results in this view.",style="font-weight:bold;font-size:24px;color:#303133 !important;font-family: OpenSans-Semibold, Helvetica Neue, Helvetica, Arial, sans-serif !important;"),
                              uiOutput(ns("conditional_scroll_box"))),
                          
                          wellPanel(
                            div(colourpicker::colourInput(ns("col1"),"Select color","#1B9E77"),
                                colourpicker::colourInput(ns("col2"),"Select color","#D95F02"),
                                colourpicker::colourInput(ns("col3"),"Select color","#7570B3"),
                                colourpicker::colourInput(ns("col4"),"Select color","#E7298A"),
                                colourpicker::colourInput(ns("col5"),"Select color","#66A61E"),
                                # colourpicker::colourInput("col6","Select color","#E6AB02"),
                                # colourpicker::colourInput("col7","Select color","#A6761D"),
                                fluidRow(column(6,textInput(ns("plot1width"),"Plot width",value = 400)),
                                         column(6,textInput(ns("plot1height"),"Plot height",value = 400))))),
                          cellWidths = c("70%","30%"))
              
    )
  )
}



overview_mainServer <- function(id) {
  moduleServer(id, function(input, output,session) {
    data_list<-reactiveValues()
    observeEvent(input$over_progress,{
      withProgress(message = 'processing...', value = 0, {
        infile <- input$file
        if (is.null(infile)) {
          return(NULL)
        }
        load(infile$datapath)
        df.pheno <- marvel$SplicePheno
        cell_meta_type<-as.data.frame(table(df.pheno$cell.type))
        data<-list()
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
          # Tabulate expressed events
          marvel <- CountEvents(MarvelObject = marvel,
                                sample.ids = sample.ids,
                                min.cells = num
          )
          data[[j]]<-marvel$N.Events$Table
          names(data)[j]<-i
          j<-j+1
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
      # browser()
      data_list$over=data},
      ignoreInit = TRUE)
    
    
    
    radarplot<-reactive({
      req(length(data_list$over) > 0)
      chart_outputs <- lapply(seq_along(data_list$over), function(i) {
        chart_id <- paste("chart", i, sep = "_")
        withSpinner(highchartOutput(session$ns(chart_id), height = '400px'),type=8)
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
    
    
    observeEvent(input$over_action,{
      warm_colors <- c(input$col1,input$col2,input$col3,input$col4,input$col5,input$col6,input$col7)
      lapply(seq_along(data_list$over), function(i) {
        chart_id <- paste("chart", i, sep = "_")
        data <- data_list$over[[i]]
        output[[chart_id]] <- renderHighchart({
          
          
          highchart() %>%
            hc_chart(backgroundColor = "#FFFFFF",type = "pie") %>%
            hc_plotOptions(pie = list(innerSize = '60%')) %>%
            hc_colors(warm_colors) %>% 
            hc_series(list(
              name = 'freq',
              data = data %>% mutate(name = event_type, y =freq) %>% select(name, y) %>% list_parse()
            )) %>%
            hc_exporting(enabled = TRUE,sourceWidth = input$plot1width, sourceHeight = input$plot1height) %>%
            
            # hc_size(height = 600, width = 800)%>%
            hc_title(text = names(data_list$over)[i])
        })
      })
      
      output$conditional_scroll_box <- renderUI({
        if (length(data_list$over) > 0) {
          div(style ="overflow-y: scroll; max-height: 800px; border: 1px solid #ddd; padding: 20px;",
              uiOutput(session$ns("charts"))
          )
        }
      })
      output$charts<-renderUI({radarplot()})
    })
    
  }
  )
}