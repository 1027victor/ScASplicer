Differential_UI<-function(id){
  ns<-NS(id)

sidebarLayout(
  sidebarPanel(
    
    
    fileInput(ns("DA_file"), label = p("File input:",style="font-weight: bold; color: #000000; font-size: 16px;")),
    wellPanel(p("Differential gene expression analysis parameter",style="font-weight: bold; color: #000000; font-size: 16px;"),
              sliderTextInput(ns("pvalue_slider"),
                              label = 'Threshold of P.Value.adj',
                              grid = TRUE,
                              choices = c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.3, 0.5, 1),
                              selected = 0.1),
              numericInput(ns("log2fc_threshold"), 'log2FC', 3,1,10),
              numericInput(ns("log2fc_threshold_"), '-log2FC', -3,-10,-1)),
    
    br(),
    wellPanel(p("Differential (spliced) gene analysis parameter",style="font-weight: bold; color: #000000; font-size: 16px;"),
              fluidRow(column(6,sliderTextInput(ns("pvalue_ad_slider"),
                              label = 'Threshold of P.Value.adj',
                              grid = TRUE,
                              choices = c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.3, 0.5, 1),
                              selected = 0.1)),
                       column(6,sliderTextInput(ns("pvalue_dts_slider"),
                              label = 'Threshold of P.Value.adj',
                              grid = TRUE,
                              choices = c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.3, 0.5, 1),
                              selected = 0.1))),
                      p("Subset sig events",style="font-weight: bold; color: #000000; font-size: 14px;"),       
                      numericInput(ns("psi_delta"), 'psi.delta', 0,1,100),
                      numericInput(ns("log2fc_psi_threshold"), 'log2FC', 0.5,1,10),
                      numericInput(ns("log2fc_psi_threshold_"), '-log2FC', -0.5,-10,-1),
                      p("Subset sig genes",style="font-weight: bold; color: #000000; font-size: 14px;"),
                      numericInput(ns("log2fc_gene_threshold"), 'log2FC', 3,1,10),
                      numericInput(ns("log2fc_gene_threshold_"), '-log2FC', -3,-10,-1),
                      numericInput(ns("logp_val_adj_threshold"),'-log10(p.val.adj)',5,1,100)),
    br(),
    fluidRow( column(4,actionButton(ns("progress"), "calculate",class='btn-info')),
              column(4,actionButton(ns("action"), "Visualization",class = "btn-primary",icon=icon('square-poll-horizontal'))))
    
  ),
  mainPanel(splitLayout(div(
    tabsetPanel(
      
      
      tabPanel("Differential gene expression analysis",icon=icon("arrow-alt-circle-up", lib = "font-awesome"),
               helpText("Show the results in this view."),
               uiOutput(ns("DA_gene_conditional_scroll_box")),
               
               helpText("Show the output_data in this view."),
               
               uiOutput(ns("DA_ex2")),
               
               DT::dataTableOutput(ns('ex2'))
               ),
      
      
      tabPanel("Differential splice analysis",icon=icon("scissors", lib = "font-awesome"),
               helpText("Show the results in this view."), uiOutput(ns("DA_splice_conditional_scroll_box")),helpText("Show the output_data in this view."),
               uiOutput(ns("DA_splice_ex3")),
               DT::dataTableOutput(ns('ex3'))
      ),
      
      tabPanel("Differential (spliced) gene analysis", icon=icon("chart-line", lib = "font-awesome"),
               helpText("Show the results in this view."), uiOutput(ns("DA_spliced_gene_conditional_scroll_box")),helpText("Show the output_data in this view."),
               uiOutput(ns("DA_spliced_gene_ex4")),
               DT::dataTableOutput(ns('ex4'))
      )   
      
    )),  wellPanel(div(
      # sliderTextInput(ns("pvalue_slider"),
      #                                  label = 'Threshold of P.Value.adj',
      #                                  grid = TRUE,
      #                                  choices = c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.3, 0.5, 1),
      #                                  selected = 0.1),
                       # numericInput(ns("log2fc_threshold"), 'log2FC', 3,1,10),
                       # numericInput(ns("log2fc_threshold_"), '-log2FC', -3,-10,-1),
                       # numericInput(ns("logp_val_adj_threshold"),'-log10(p.val.adj)',5,1,100),
                       numericInput(ns("radius_num"),'radius',1,1,10),                                  
                       fluidRow(
                         column(4, colourpicker::colourInput(ns("colorUp"), 'Up', "#450456")),
                         column(4, colourpicker::colourInput(ns("colorStable"), 'Stable', "#FDE93A")),
                         column(4, colourpicker::colourInput(ns("colorDown"),'Down', "#59ACA9"))),
                       
                       p("Differential gene expression analysis plot parameter",style="font-weight: bold; color: #000000; font-size: 16px;white-space: pre-wrap; word-wrap: break-word;"),
                       fluidRow(column(6,textInput(ns("plot1width_DA"),"Plot width",value = 400)),
                                column(6,textInput(ns("plot1height_DA"),"Plot height",value = 400))),
                      p("Differential gene expression analysis plot parameter",style="font-weight: bold; color: #000000; font-size: 16px;white-space: pre-wrap; word-wrap: break-word;"),      
                      fluidRow(column(6,textInput(ns("plot1width_DS"),"Plot width",value = 1500)),
                              column(6,textInput(ns("plot1height_DS"),"Plot height",value = 300))),        
                      p("Differential gene expression analysis plot parameter",style="font-weight: bold; color: #000000; font-size: 16px;white-space: pre-wrap; word-wrap: break-word;"),
                      fluidRow(column(6,textInput(ns("plot1width_DSG"),"Plot width",value = 400)),
                                column(6,textInput(ns("plot1height_DSG"),"Plot height",value = 400)))
                                
                                )),
    
    cellWidths = c("75%","25%"))
  )
)
}



Differential_mainServer <- function(id) {
  moduleServer(id, function(input, output,session) {

data_list_DA <- reactiveValues()

# process data
observeEvent(input$progress,{
  # print("Button clicked")
  withProgress(message = 'processing...', value = 0, {
    infile <- input$DA_file
    if (is.null(infile)) {
      return(NULL)
    }
    load(infile$datapath)
    df.pheno <- marvel$SplicePheno
    data<-list()
    data_splice<-list()
    data_splice_genes<-list()
    # merge_df<-NULL
    merge_df_Global<-NULL
    merge_df_splice<-NULL
    merge_df_spliced_genes<-NULL
    l<-1
    cell_types<-unique(df.pheno$cell.type)
    cell_meta_type<-as.data.frame(table(df.pheno$cell.type))
    # min_cell
    num<-round(min(cell_meta_type$Freq/2))+1
    if(num>25){num=25}
    total_iterations <- (length(cell_types) - 1) * length(cell_types) / 2
    progress <- 0
    # browser()
    # print(total_iterations)
    # withProgress(message = 'processing...', value = 0, {
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
        marvel <- PlotDEValues(MarvelObject=marvel,
                               pval=input$pvalue_slider,
                               log2fc=0.5,
                               point.size=0.1,
                               xlabel.size=8,
                               level="gene.global",
                               anno=FALSE
        )
      
        # Differential splicing analysis
        # num<-min(num_one,num_two)
        # num<-round(min(num_one,num_two)/4)+1
        # # print(num)
        # if(num>25){num=25}
        # else{num<-round(min(num_one,num_two)/4)+1}
        marvel <- CompareValues(MarvelObject=marvel,
                                cell.group.g1=cell.group.g1,
                                cell.group.g2=cell.group.g2,
                                min.cells=num,
                                method=c("ad", "dts"),
                                method.adjust="fdr",
                                level="splicing",
                                assign.modality = FALSE,
                                event.type=c("SE", "MXE", "RI", "A5SS", "A3SS"),
                                show.progress=FALSE
        )
        # browser()
        
        
        #Differential (spliced) gene analysis
        marvel <- CompareValues(MarvelObject=marvel,
                                cell.group.g1=cell.group.g1,
                                cell.group.g2=cell.group.g2,
                                psi.method=c("ad", "dts"),
                                #set shiny parameter
                                psi.pval=c(input$pvalue_ad_slider, input$pvalue_dts_slider),
                                #set shiny parameter
                                psi.delta=input$psi_delta,
                                method.de.gene="wilcox",
                                method.adjust.de.gene="fdr",
                                downsample=FALSE,
                                show.progress=FALSE,
                                level="gene.spliced"
        )
        df_spliced_genes<-marvel$DE$Exp.Spliced$Table
        # df<-marvel$DE$Exp$Table
        df_global <- marvel$DE$Exp.Global$Table
        # data_splice_event<-marvel$DE$PSI$Table
        df_ad<-marvel$DE$PSI$Table[["ad"]]
        df_dts<-marvel$DE$PSI$Table[["dts"]]
        data_splice_event<-rbind(df_ad,df_dts)
        data_splice_event<-data_splice_event[,c("tran_id","gene_id", "gene_short_name", "event_type","gene_type","mean.diff","mean.g1","mean.g2")]
        data_splice_event<-unique(data_splice_event)
        #set shiny parameter
        index_ad <- which(abs(df_ad$mean.diff) > input$psi_delta & df_ad$p.val.adj < input$pvalue_ad_slider & df_ad$outliers==FALSE)
        index_dts <- which(abs(df_dts$mean.diff) > input$psi_delta & df_dts$p.val.adj < input$pvalue_dts_slider & df_dts$outliers==FALSE)
        
        df_ad<-df_ad[index_ad,]
        df_dts<-df_dts[index_dts,]
        df_splice<-rbind(df_ad,df_dts)
        
        df_splice_event_genes<-df_splice[,c("gene_id", "gene_short_name", "gene_type")]
        df_splice_event_genes<-unique(df_splice_event_genes)
        
        # Annotate gene pval, log2fc
        df_spliced_genes_plot <- left_join(df_splice_event_genes, df_spliced_genes[,c("gene_id", "log2fc", "p.val.adj")], by="gene_id")
        df_spliced_genes_plot$log2fc[is.na(df_spliced_genes_plot$log2fc)] <- 0
        df_spliced_genes_plot$p.val.adj[is.na(df_spliced_genes_plot$p.val.adj)] <- 1
        
        # Indicate sig events and direction
        #set shiny parameter
        df_spliced_genes_plot$sig <- NA
        df_spliced_genes_plot$sig[which(df_spliced_genes_plot$p.val.adj < input$pvalue_ad_slider&df_spliced_genes_plot$log2fc > input$log2fc_psi_threshold)] <- "up"
        df_spliced_genes_plot$sig[which(df_spliced_genes_plot$p.val.adj < input$pvalue_ad_slider&df_spliced_genes_plot$log2fc < input$log2fc_psi_threshold_)] <- "down"
        df_spliced_genes_plot$sig[is.na(df_spliced_genes_plot$sig)] <- "n.s."
        # df$sig <- factor(df$sig, levels=c("up", "down", "n.s."))
        
        
        
        df_splice<-df_splice[,c('event_type','mean.g1','mean.g2','gene_short_name')]
        df_splice<-unique(df_splice)
        df_splice$gene_short_name.event_type <- paste(df_splice$gene_short_name,"_",df_splice$event_type,sep="")
        df_splice<-df_splice[,c('gene_short_name.event_type','mean.g1','mean.g2','gene_short_name')]
        
        # browser()
        # trans long data
        data_long <- melt(df_splice[,c('gene_short_name.event_type','mean.g1','mean.g2')], id.vars = "gene_short_name.event_type")
        merged_data <- left_join(data_long, select(df_splice,gene_short_name.event_type, gene_short_name), by = "gene_short_name.event_type")
        
        # browser()
        df_spliced_genes_plot$group_one<-paste(cell_types[i])
        df_spliced_genes_plot$group_two<-paste(cell_types[j])
        data_splice_event$group_one<-paste(cell_types[i])
        data_splice_event$group_two<-paste(cell_types[j])
        # df$comparison<-paste(cell_types[i], " vs ", cell_types[j])
        # df_global$comparsion<-paste(cell_types[i], " vs ", cell_types[j])
        df_global$group_one<-paste(cell_types[i])
        df_global$group_two<-paste(cell_types[j])
        
        
        # DA data
        data[[l]]<-marvel$DE$Exp.Global$Table
        data[[l]]$z <- as.factor(data[[l]]$sig)
        data[[l]]$label <- ifelse(data[[l]]$log2fc > input$log2fc_threshold | data[[l]]$log2fc < input$log2fc_threshold_, data[[l]]$gene_short_name, "")

        data_splice[[l]]<-merged_data
        # browser()
        data_splice_genes[[l]]<-df_spliced_genes_plot
        data_splice_genes[[l]]$z <- as.factor(data_splice_genes[[l]]$sig)
        data_splice_genes[[l]]$label <- ifelse((data_splice_genes[[l]]$log2fc > input$log2fc_gene_threshold |data_splice_genes[[l]]$log2fc < input$log2fc_gene_threshold_) & -log10(data_splice_genes[[l]]$p.val.adj)>input$logp_val_adj_threshold, data_splice_genes[[l]]$gene_short_name, "")
    

        names(data_splice)[l]<-paste(cell_types[i], " vs ", cell_types[j])
        # browser()
        names(data)[l]<-paste(cell_types[i], " vs ", cell_types[j])
        names(data_splice_genes)[l]<-paste(cell_types[i], " vs ", cell_types[j])
        l<-l+1
        if(is.null(merge_df_Global)) {
          merge_df_Global <- df_global
        } else {
          merge_df_Global <- rbind(merge_df_Global, df_global)
        }
        # if(is.null(merge_df)) {
        #   merge_df <- df
        # } else {
        #   merge_df <- rbind(merge_df, df)
        # }
        if(is.null(merge_df_splice)) {
          merge_df_splice <- data_splice_event
        } else {
          merge_df_splice <- rbind(merge_df_splice, data_splice_event)
        }
        
        if(is.null(merge_df_spliced_genes)) {
          merge_df_spliced_genes <- df_spliced_genes_plot
        } else {
          merge_df_spliced_genes <- rbind(merge_df_spliced_genes, df_spliced_genes_plot)
        }
        
        progress <- progress + 1
        setProgress(value = progress / total_iterations,message = sprintf("Processing progress: %d%%", round((progress / total_iterations) * 100)))
      }
    }
    # showNotification("complete! ",closeButton = TRUE,
    #                  duration = 60,type = c("message"))
    # showModal(modalDialog(
    #   "complete!",
    #   easyClose = TRUE,
    #   footer = modalButton("close")
    # ))
    
    showModal(modalDialog(
      title = "",
      "The data processing is complete!",
      footer = modalButton("close")
    ))
    
  })
  # browser()
  data_list_DA$data_DA=data
  # data_list_DA$data_merge=merge_df
  data_list_DA$data_merge_global=merge_df_Global
  data_list_DA$data_merge_splice=merge_df_splice
  data_list_DA$data_spliced=data_splice
  data_list_DA$data_merge_spliced_genes= merge_df_spliced_genes
  data_list_DA$data_spliced_genes=data_splice_genes
  
  # list(data_DA=data,data_merge=merge_df,data_merge_global=merge_df_Global)
  # print(data_DA)
})


radarplot_DA<-reactive({
  req(length(data_list_DA$data_DA) > 0)
  chart_outputs <- lapply(seq_along(data_list_DA$data_DA), function(i) {
    chart_id <- paste("chart_DA", i, sep = "_")
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



# splice chart
radarplot_DA_spliced<-reactive({
  req(length(data_list_DA$data_spliced) > 0)
  chart_outputs <- lapply(seq_along(data_list_DA$data_spliced), function(i) {
    data<-data_list_DA$data_spliced[[i]]
    base_height_per_row <- 120
    base_width_per_col <- 40
    rows <- 2
    cols <-nrow(data)
    height <- base_height_per_row * rows
    width <- base_width_per_col * cols
    chart_id <- paste("chart_spliced", i, sep = "_")
    withSpinner(highchartOutput(session$ns(chart_id),width=paste(width,"px",sep=''),height = paste(height,"px",sep='')),type=8)
  })
  
  chart_rows <- lapply(seq(1, length(chart_outputs), by = 1), function(i) {
    fluidRow(
      column(12,chart_outputs[[i]])
      # if (i + 1 <= length(chart_outputs)) column(6, chart_outputs[[i + 1]])
      # if (i + 1 <= length(chart_outputs)) column(4, chart_outputs[[i + 2]])
    )
  })
  do.call(tagList, chart_rows)
})



#spliced_genes
radarplot_spliced_genes<-reactive({
  req(length(data_list_DA$data_spliced_genes) > 0)
  chart_outputs <- lapply(seq_along(data_list_DA$data_spliced_genes), function(i) {
    chart_id <- paste("chart_spliced_genes", i, sep = "_")
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






values <- reactiveValues()

# gene plot
observeEvent(input$action,{
  
  values$plot1height_DA<-input$plot1height_DA
  values$plot1width_DA<-input$plot1width_DA
  values$plot1height_DS<-input$plot1height_DS
  values$plot1width_DS<-input$plot1width_DS
  values$plot1height_DSG<-input$plot1height_DSG
  values$plot1width_DSG<-input$plot1width_DSG
  values$colorUp<-input$colorUp
  values$colorDown<-input$colorDown
  values$colorStable<-input$colorStable
  values$radius_num<-input$radius_num

  lapply(seq_along(data_list_DA$data_DA), function(i) {
    data <- data_list_DA$data_DA[[i]]
    # data$z <- as.factor(data$sig)
    # data$label <- ifelse(data$log2fc > input$log2fc_threshold | data$log2fc < input$log2fc_threshold_, data$gene_short_name, "")
    # setting plot color
    color_mapping <- setNames(c(values$colorUp,values$colorDown,values$colorStable), levels(data$z))
    
    data$color <- color_mapping[data$z]
    chart_id <- paste("chart_DA", i, sep = "_")
    output[[chart_id]] <- renderHighchart({
      
      hchart(data, "scatter", hcaes(x = log2fc, y = -log10(p.val.adj), color = color, label = label)) %>%
        hc_plotOptions(
          scatter = list(
            boostThreshold = 5000,
            animation = FALSE,
            dataLabels = list(
              enabled = TRUE,
              format = '{point.label}',
              style = list(textOutline = FALSE),
              allowOverlap = FALSE
              
            ),
            marker = list(
              radius = values$radius_num
            ),
            states = list(
              hover = list(
                enabled = FALSE
              )
            )
          )
        ) %>%
        # hc_colors(custom_colors) %>%
        hc_xAxis(
          title = list(text = 'log2FC', style = list(fontSize = "14px", color = "black", fontWeight = "bold")),
          gridLineWidth = 0,
          lineWidth = 1,
          lineColor = "black",
          labels = list(style = list(fontSize = "14px", color = "black", fontWeight = "bold")),
          tickInterval = 2
        ) %>%
        hc_yAxis(
          title = list(text = '-log10(p-value)', style = list(fontSize = "14px", color = "black", fontWeight = "bold")),
          gridLineWidth = 0,
          lineWidth = 1,
          lineColor = "black",
          labels = list(style = list(fontSize = "14px", color = "black", fontWeight = "bold")),
          tickInterval = 5
        ) %>%
        hc_title(text = names(data_list_DA$data_DA)[i]) %>%
        hc_legend(enabled = FALSE) %>%
        # hc_tooltip(useHTML = TRUE, headerFormat = "", pointFormat = '{point.gene_short_name}') %>%
        hc_tooltip(enabled=FALSE) %>%
        hc_chart(
          backgroundColor = "#FFFFFF",
          style = list(fontFamily = 'Arial')
        ) %>%
        hc_exporting(enabled = TRUE,sourceWidth = values$plot1width_DA, sourceHeight = values$plot1height_DA)
    })
    
  })
  
  output$DA_gene_conditional_scroll_box <- renderUI({
    if (length(data_list_DA$data_DA) <0){
      return (NULL)}
    else
    {div(style ="overflow-y: scroll; max-height: 800px;padding: 20px;",
         uiOutput(session$ns("charts_DA"))
    )}
  })
  
  output$charts_DA<-renderUI({radarplot_DA()})
  
  
  output$DA_ex2 <- renderUI({
    if (is.null(data_list_DA)) {
      return(NULL)
    }
    div(class = "dropdown",
        tags$button(class = "btn btn-custom dropdown-toggle", type = "button", id = "dropdownMenuButton",
                    `data-bs-toggle` = "dropdown", `aria-haspopup` = "true", `aria-expanded` = "false",
                    "Download  Data ", span(class = "caret")
        ),
        div(class = "dropdown-menu", `aria-labelledby` = "dropdownMenuButton",
            downloadLink(session$ns("download_csv"), "CSV", class = "btn-link"),
            downloadLink(session$ns("download_txt"), "TXT", class = "btn-link")
        ))
    
  })
  
  
  output$ex2 <- DT::renderDataTable({
    if (is.null(data_list_DA)){
      return(NULL)}
    DT::datatable(data_list_DA$data_merge_global,
                  options = list(pageLength =10, scrollX = TRUE))}
  )
  
  output$download_csv <- downloadHandler(
    filename = function() {
      paste("result", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      # load your_data_frame 
      write.csv(data_list_DA$data_merge_global, file)
    }
  )
  
  # downlaoad txt
  output$download_txt <- downloadHandler(
    filename = function() {
      paste("result", Sys.Date(), ".txt", sep="")
    },
    content = function(file) {
      # load your_data_frame
      write.table(data_list_DA$data_merge_global, file, row.names = FALSE)
    }
  )
  
  
  #splice event plot
  lapply(seq_along(data_list_DA$data_spliced), function(i) {
    data <- data_list_DA$data_spliced[[i]]
    chart_id <- paste("chart_spliced", i, sep = "_")
    output[[chart_id]] <- renderHighchart({
      
      base_height_per_row <- 120
      base_width_per_col <- 40
      rows <- 2
      cols <- nrow(data)
      height <- base_height_per_row * rows
      width <- base_width_per_col * cols
      # hc_chart(
      #   backgroundColor = "#FFFFFF",
      #   style = list(fontFamily = 'Arial')
      # ) %>%
      hchart(data, "heatmap", hcaes(x = gene_short_name.event_type, y = variable, value = value)) %>%
        hc_chart(width = width,height=height)%>%
        # hc_chart(width = 1000) %>%  # load width
        # hc_xAxis(title = list(text = NULL), labels = list(style = list(fontSize = "14px", color = "black"),rotation = 45))%>% 
        # hc_xAxis(title = list(text = ""), labels = list(style = list(fontSize = "14px", color = "black")))%>%
        hc_xAxis(
          title = list(text = ""),
          labels = list(style = list(fontSize = "12px", color = "black")),
          rotation = -45
        )%>%
        hc_yAxis(title = list(text = ""), labels = list(style = list(fontSize = "12px", color = "black")),
                 categories = c((strsplit(names(data_list_DA$data_spliced)[i],split=" ")[[1]])[1], (strsplit(names(data_list_DA$data_spliced)[i],split=" ")[[1]])[5]) 
                 )%>%
        
        hc_tooltip(useHTML = TRUE, headerFormat = "", pointFormat = '{point.gene_short_name}:{point.value:.2f}') %>%
        # hc_chart(width = 1600, height = 400)%>%
        hc_legend(title = list(text = 'Mean PSI'), align = 'center', layout = 'horizontal', verticalAlign = 'top') %>%
        hc_colorAxis(stops = color_stops(n =nrow(data) , colors = c("#3060cf", "#fffbbc", "#c4463a")))%>%
        # hc_chart(backgroundColor = "#FFFFFF",style = list(fontFamily = 'Arial'))%>%
        hc_chart(events = list(load = JS("function() {
            var chart = this;
            chart.update({
              chart: {
                backgroundColor: '#FFFFFF' // set explort background color
              }
            });
            console.log('Updated chart background color!');
          }")))%>%
        hc_exporting(enabled = TRUE,sourceWidth = values$plot1width_DS, sourceHeight = values$plot1height_DS)
      # hc_chart(backgroundColor = "#FFFFFF",style = list(fontFamily = 'Arial')) 
    })
    
  })
  
  output$DA_splice_conditional_scroll_box <- renderUI({
    if (length(data_list_DA$data_spliced) <0){
      return (NULL)}
    else
    {div(style ="overflow-y: scroll;overflow-x:scroll;max-height: 800px;padding: 20px;",
         # highchartOutput(session$ns("charts_DA_spliced"),width='1000px',height = '300px')
         radarplot_DA_spliced()
    )}
  })
  
  # output$charts_DA_spliced<-renderUI({radarplot_DA_spliced()})
  
  output$DA_splice_ex3 <- renderUI({
    if (is.null(data_list_DA)) {
      return(NULL)
    }
    div(class = "dropdown",
        tags$button(class = "btn btn-custom dropdown-toggle", type = "button", id = "dropdownMenuButton",
                    `data-bs-toggle` = "dropdown", `aria-haspopup` = "true", `aria-expanded` = "false",
                    "Download Data ", span(class = "caret")
        ),
        div(class = "dropdown-menu", `aria-labelledby` = "dropdownMenuButton",
            downloadLink(session$ns("download_PSI_csv"), "CSV", class = "btn-link"),
            downloadLink(session$ns("download_PSI_txt"), "TXT", class = "btn-link")
        ))
    
  })
  
  
  output$ex3 <- DT::renderDataTable({
    if (is.null(data_list_DA)){
      return(NULL)}
    DT::datatable(data_list_DA$data_merge_splice,
                  options = list(pageLength =10, scrollX = TRUE))}
  )
  
  output$download_PSI_csv <- downloadHandler(
    filename = function() {
      paste("result", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      # load your_data_frame 
      write.csv(data_list_DA$data_merge_splice, file)
    }
  )
  
  # downlaoad txt
  output$download_PSI_txt <- downloadHandler(
    filename = function() {
      paste("result", Sys.Date(), ".txt", sep="")
    },
    content = function(file) {
      # load your_data_frame 
      write.table(data_list_DA$data_merge_splice, file, row.names = FALSE)
    }
  )
  
  
  #spliced_genes
  lapply(seq_along(data_list_DA$data_spliced_genes), function(i) {
    data <- data_list_DA$data_spliced_genes[[i]]
    # data$z <- as.factor(data$sig)
    # data$label <- ifelse((data$log2fc > input$log2fc_threshold | data$log2fc < input$log2fc_threshold_) & -log10(data$p.val.adj)>input$logp_val_adj_threshold,data$gene_short_name, "")
    
    color_mapping <- setNames(c(values$colorUp,values$colorDown,values$colorStable), levels(data$z))
    
    
    data$color <- color_mapping[data$z]
    chart_id <- paste("chart_spliced_genes", i, sep = "_")
    output[[chart_id]] <- renderHighchart({
      hchart(data, "scatter", hcaes(x = log2fc, y = -log10(p.val.adj), color = color, label = label)) %>%
        hc_plotOptions(
          scatter = list(
            boostThreshold = 5000,
            animation = FALSE,
            dataLabels = list(
              enabled = TRUE,
              format = '{point.label}',
              style = list(textOutline = FALSE),
              allowOverlap = FALSE
              
            ),
            marker = list(
              radius = values$radius_num
            ),
            states = list(
              hover = list(
                enabled = FALSE
              )
            )
          )
        ) %>%
        # hc_colors(custom_colors) %>%
        hc_xAxis(
          title = list(text = 'log2FC', style = list(fontSize = "14px", color = "black", fontWeight = "bold")),
          gridLineWidth = 0,
          lineWidth = 1,
          lineColor = "black",
          labels = list(style = list(fontSize = "14px", color = "black", fontWeight = "bold")),
          tickInterval = 2
        ) %>%
        hc_yAxis(
          title = list(text = '-log10(p-value)', style = list(fontSize = "14px", color = "black", fontWeight = "bold")),
          gridLineWidth = 0,
          lineWidth = 1,
          lineColor = "black",
          labels = list(style = list(fontSize = "14px", color = "black", fontWeight = "bold")),
          tickInterval = 5
        ) %>%
        hc_title(text = names(data_list_DA$data_DA)[i]) %>%
        hc_legend(enabled = FALSE) %>%
        # hc_tooltip(useHTML = TRUE, headerFormat = "", pointFormat = '{point.gene_short_name}') %>%
        hc_tooltip(enabled=FALSE) %>%
        hc_chart(
          backgroundColor = "#FFFFFF",
          style = list(fontFamily = 'Arial')
        ) %>%
        hc_exporting(enabled = TRUE,sourceWidth = values$plot1width_DSG, sourceHeight = values$plot1height_DSG)
    })
    
  })
  
  output$DA_spliced_gene_conditional_scroll_box <- renderUI({
    if (length(data_list_DA$data_spliced_genes) <0){
      return (NULL)}
    else
    {div(style ="overflow-y: scroll; max-height: 800px;padding: 20px;",
         uiOutput(session$ns("charts_spliced_genes"))
    )}
  })
  
  output$charts_spliced_genes<-renderUI({radarplot_spliced_genes()})
  
  
  output$DA_spliced_gene_ex4 <- renderUI({
    if (is.null(data_list_DA)) {
      return(NULL)
    }
    div(class = "dropdown",
        tags$button(class = "btn btn-custom dropdown-toggle", type = "button", id = "dropdownMenuButton",
                    `data-bs-toggle` = "dropdown", `aria-haspopup` = "true", `aria-expanded` = "false",
                    "Download  Data ", span(class = "caret")
        ),
        div(class = "dropdown-menu", `aria-labelledby` = "dropdownMenuButton",
            downloadLink(session$ns("download_spliced_genes_csv"), "CSV", class = "btn-link"),
            downloadLink(session$ns("download_spliced_genes_txt"), "TXT", class = "btn-link")
        ))
    
  })
  
  
  output$ex4 <- DT::renderDataTable({
    if (is.null(data_list_DA)){
      return(NULL)}
    DT::datatable(data_list_DA$data_merge_spliced_genes,
                  options = list(pageLength =10, scrollX = TRUE))}
  )
  
  output$download_spliced_genes_csv <- downloadHandler(
    filename = function() {
      paste("result", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      # load your_data_frame 
      write.csv(data_list_DA$data_merge_spliced_genes, file)
    }
  )
  
  # downlaoad txt
  output$download_spliced_genes_txt <- downloadHandler(
    filename = function() {
      paste("result", Sys.Date(), ".txt", sep="")
    },
    content = function(file) {
      # load your_data_frame
      write.table(data_list_DA$data_merge_spliced_genes, file, row.names = FALSE)
    }
  )
  
  
})

}

)

}