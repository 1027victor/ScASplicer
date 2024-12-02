my_bar<-function(data,
                 color_min = "#e1e8f0",
                 color_max = "#326fa4",
                 width=850,
                 height=600,
                 go_term=10,
                 fontSize = "12px",
                 sourceWidth = 850,
                 sourceHeight = 600){
  color_min <- color_min
  color_max <- color_max
  data=data[1:go_term,]
  data$color <- scales::col_numeric(palette = c(color_min, color_max),
                                         domain = range(data$neg_log10pvalue))(data$neg_log10pvalue)
  # browser()
  barplot<-hchart(data, "bar", hcaes(x = Description, y = Count, color = color)) %>%
    hc_chart(width = width,height=height) %>%
    # hc_title(text = "GO Terms Enrichment Barplot", align = "center") %>%
    hc_xAxis(title = list(text = " "),
             lineColor = "black",
             lineWidth = 1,
             gridLineWidth = 0,
             labels = list(
               style = list(fontSize = fontSize, color = "black", fontWeight = "normal")
             )) %>%
    hc_yAxis(title = list(text = "Gene Number",style = list(color = "black", fontSize = fontSize,fontWeight = "normal")),
             lineColor = "black",
             lineWidth = 1,
             gridLineWidth = 0,
             tickInterval = 5,
             labels = list(
               style = list(fontSize = fontSize, color = "black", fontWeight = "normal")
             )) %>%
    hc_colorAxis(stops = color_stops(length(data$neg_log10pvalue), c(color_min,color_max))) %>%
    hc_legend(title = list(text = "-log10(pvalue)",style = list(color = "black", fontSize = fontSize,fontWeight = "normal")), layout = "vertical", align = "right", verticalAlign = "middle") %>%
    hc_tooltip(pointFormat = "<b>{point.y}</b> genes") %>%
    hc_chart(events = list(load = JS("function() {
            var chart = this;
            chart.update({
              chart: {
                backgroundColor: '#FFFFFF'
              }
            });
            console.log('Updated chart background color!');
          }")))%>%
    hc_exporting(enabled = TRUE,sourceWidth = sourceWidth,sourceHeight = sourceHeight)%>%
  hc_add_theme(hc_theme_smpl())
  return(barplot)
}
# my_Bubble_plot<-function(data,
#                         width = 850,
#                         height = 600,
#                         go_term = 10,
#                         sourceWidth = 850,
#                         sourceHeight = 600,
#                         color_min = "#FFF6A5", 
#                         color_max = "#c1021f",
#                         fontSize = "12px",
#                         minsize = 10,
#                         maxsize = 30)
#                         {
#   color_min <- color_min 
#   color_max <- color_max 
#   data=data[1:go_term,]
#   data$color <- scales::col_numeric(palette = c(color_min, color_max),
#                                        domain = range(data$neg_log10pvalue))(data$neg_log10pvalue)
#   data$category_index <- 0:(length(data$Description)-1)
#   bubble_chart<-highchart() %>%
#     hc_chart(width = width,height=height) %>% 
#     hc_chart(type = 'bubble') %>%
#     hc_title(text = " ") %>%
#     hc_xAxis(title = list(text = "Gene Number",style = list(color = "black", fontSize = fontSize,fontWeight = "normal")),
#              lineColor = "black",
#              lineWidth = 1,
#              # tickInterval = 5, 
#              gridLineWidth = 1,
#              labels = list(
#                style = list(fontSize = fontSize, color = "black", fontWeight = "normal"))) %>%
#     hc_yAxis(title = list(text = ""),
#              categories = data$Description,
#              min = 0,max=length(data$Description)-1,
#              lineColor = "black",
#              lineWidth = 1,
#              gridLineWidth = 1,
#              labels = list(
#                style = list(fontSize = fontSize, color = "black", fontWeight = "normal"))) %>%
#     hc_add_series(
#       data = data,
#       type = "bubble",
#       hcaes(x = Count, y = category_index,z=sqrt(Count),color=color),
#       showInLegend = TRUE,
#       name='gene number'
#     ) %>%
#     # minSize = 10, maxSize = 30 set the size of the bubble
#     hc_plotOptions(bubble = list(minSize = minsize, maxSize = maxsize, fillOpacity = 0.2)) %>%
#     hc_colorAxis(stops = color_stops(length(data$neg_log10pvalue), c(color_min,color_max))) %>%
#     hc_legend(title = list(text = "-log10(pvalue)",style = list(color = "black", fontSize = fontSize,fontWeight = "normal")), layout = "vertical", align = "right", verticalAlign = "middle") %>%
#     hc_tooltip(pointFormat = "<b>{point.x}</b> genes") %>%
#     hc_chart(events = list(load = JS("function() {
#               var chart = this;
#               chart.update({
#                 chart: {
#                   backgroundColor: '#FFFFFF' // set the background color of the chart
#                 }
#               });
#               console.log('Updated chart background color!')
#             }")))%>%
#     hc_exporting(enabled = TRUE,sourceHeight=sourceHeight,sourceWidth=sourceWidth)
#     return(bubble_chart)
  
# }

my_Bubble_plot <- function(data,
                           width = 850,
                           height = 600,
                           go_term = 10,
                           sourceWidth = 850,
                           sourceHeight = 600,
                           color_min = "#FFF6A5",
                           color_max = "#c1021f",
                           fontSize = "12px",
                           minsize = 10,
                           maxsize = 30) {
  color_min <- color_min
  color_max <- color_max
  data = data[1:go_term,]
  # 计算GeneRatio
  data$GeneRatio <- sapply(strsplit(as.character(data$GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  # 按GeneRatio排序
  data = data %>% arrange(GeneRatio)
  data$category_index <- 0:(length(data$Description) - 1)
  
  bubble_chart <- highchart() %>%
    hc_chart(width = width, height = height) %>%
    hc_chart(type = 'bubble') %>%
    hc_title(text = " ") %>%
    hc_xAxis(title = list(text = "GeneRatio", style = list(color = "black", fontSize = fontSize, fontWeight = "normal")),
             lineColor = "black",
             lineWidth = 1,
             gridLineWidth = 1,
             labels = list(
               style = list(fontSize = fontSize, color = "black", fontWeight = "normal"))) %>%
    hc_yAxis(title = list(text = ""),
             categories = data$Description,
             min = 0, max = length(data$Description) - 1,
             lineColor = "black",
             lineWidth = 1,
             gridLineWidth = 1,
             labels = list(
               style = list(fontSize = fontSize, color = "black", fontWeight = "normal"))) %>%
    hc_add_series(
      data = data,
      type = "bubble",
      hcaes(x = GeneRatio, y = category_index, z = Count),  # 移除了color的映射
      showInLegend = TRUE,
      name = 'Count',
      colorKey = 'neg_log10pvalue',  # 指定颜色键
      color = 'black'
    ) %>%
    hc_plotOptions(
      bubble = list(minSize = minsize, maxSize = maxsize, fillOpacity = 0.6),  # 调整透明度以更清晰地显示颜色
      series = list(
        point = list(
          events = list(
            hide = JS("function () { this.series.chart.colorAxis[0].update({}, true); }"),
            show = JS("function () { this.series.chart.colorAxis[0].update({}, true); }")
          )
        )
      )
    ) %>%
    hc_colorAxis(
      stops = color_stops(10, c(color_min, color_max)),
      min = min(data$neg_log10pvalue, na.rm = TRUE),
      max = max(data$neg_log10pvalue, na.rm = TRUE),
      title = list(text = "-log10(pvalue)", style = list(color = "black", fontSize = fontSize, fontWeight = "normal"))
    ) %>%
    # hc_legend(layout = "vertical", align = "right", verticalAlign = "middle") %>%
    hc_legend(title = list(text = "-log10(pvalue)",style = list(color = "black", fontSize = fontSize,fontWeight = "normal")), layout = "vertical", align = "right", verticalAlign = "middle") %>%
    hc_tooltip(pointFormat = "<b>{point.z}</b> genes") %>%  # 在提示框中显示neg_log10pvalue
    hc_chart(events = list(load = JS("function() {
                  var chart = this;
                  chart.update({
                    chart: {
                      backgroundColor: '#FFFFFF' // 设置图表背景颜色
                    }
                  });
                  console.log('Updated chart background color!')
                }"))) %>%
    hc_exporting(enabled = TRUE, sourceHeight = sourceHeight, sourceWidth = sourceWidth)
  return(bubble_chart)
}

my_circle<-function(data,
                    go_term=10,
                    outcolor='#69c3c5',
                    outfontcolor='#FFFFFF',
                    upregulated='#a0c6c9',
                    downregulated='#c05678',
                    Fourthcirclecolor='#69c3c5',
                    lengend_down='#FF906F',
                    lengend_up='#861D30'
                  # Ontology='BP'
                  ){
    # data=data[1:go_term,]
    circos.clear()
    circle_size = unit(1, 'snpc')
    circos.par(gap.degree = 2, start.degree = 90)
    data=data[1:go_term,]
    plot_data <- data[c('ID', 'gene_num.min', 'gene_num.max')] 
    ko_color <- c(rep(outcolor,nrow(data))) #各二级分类的颜色和数目

    circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1) 
    #track.height默认调整轨道的宽度,#bg.border是设置扇区外框的颜色
    circos.track(
    ylim = c(0, 1), track.height = 0.1, bg.border = ko_color, bg.col = ko_color,  
    panel.fun = function(x, y) {
    ylim = get.cell.meta.data('ycenter')  
    xlim = get.cell.meta.data('xcenter')
    sector.name = get.cell.meta.data('sector.index')  
    circos.axis(h = 'top', labels.cex = 0.5, labels.niceFacing = FALSE) 
    #font是设置字体的样式
    circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE,col = outfontcolor,font=2)  
    } )


    # 第二圈，绘制富集的基因和富集p值
    plot_data <- data[c('ID', 'gene_num.min', 'gene_num.rich', 'neg_log10pvalue')]  
    label_data <- as.data.frame(data['gene_num.rich'])  
    p_max <- round(max(data$'neg_log10pvalue')) + 1  
    colorsChoice <- colorRampPalette(c(lengend_down, lengend_up))  
    color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))

    circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.12, bg.border = NA, stack = TRUE,  
    panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  
        ylim = get.cell.meta.data('ycenter')  
        xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
        sector.name = label_data[get.cell.meta.data('sector.index'),1]
        circos.text(xlim, ylim, sector.name, cex = 0.7, niceFacing = FALSE)
    } )

    # 第三圈，绘制上下调基因
    data$up <- data$Upregulated_Prop * data$gene_num.max
    plot_data_up <- data[c('ID', 'gene_num.min', 'up')]
    names(plot_data_up) <- c('ID', 'start', 'end')
    plot_data_up$type <- 1  

    data$down <- data$Downregulated_Prop * data$gene_num.max
    plot_data_down <- data[c('ID', 'gene_num.min', 'down')]
    names(plot_data_down) <- c('ID', 'start', 'end')
    plot_data_down$type <- 2  

    plot_data <- rbind(plot_data_up, plot_data_down)
    label_data <- data[c('up', 'down', 'Upregulated', 'Downregulated')]
    color_assign <- colorRamp2(breaks = c(1, 2), col = c(upregulated, downregulated))

    circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.09, bg.border = NA, stack = TRUE,
    panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = color_assign(value[[1]]), border = color_assign(value[[1]]), ...)  
        ylim = get.cell.meta.data('cell.bottom.radius') - 0.5 
        xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
        sector.name = label_data[get.cell.meta.data('sector.index'),3]
        circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE)
        xlim = (label_data[get.cell.meta.data('sector.index'),2]+label_data[get.cell.meta.data('sector.index'),1]) / 2
        sector.name = label_data[get.cell.meta.data('sector.index'),4]
        circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = FALSE)
    } )




    # 第四圈，绘制富集因子
    plot_data <- data[c('ID', 'gene_num.min', 'gene_num.max', 'rich.factor')] 
    label_data <- data['category']  
    if (unique(label_data) == 'BP') {
       color_assign <- c('BP' = Fourthcirclecolor)#各二级分类的名称和颜色
    }
    # color_assign <- c('BP' = '#69c3c5')#各二级分类的名称和颜色
    else if (unique(label_data) == 'CC') {
        color_assign <- c('CC' = Fourthcirclecolor)
    }
    else  {
        color_assign <- c('MF' = Fourthcirclecolor)
    }
    circos.genomicTrack(
    plot_data, ylim = c(0, 1), track.height = 0.3, bg.col = 'gray95', bg.border = NA,  
    panel.fun = function(region, value, ...) {
        sector.name = get.cell.meta.data('sector.index')  
        circos.genomicRect(region, value, col = color_assign[label_data[sector.name,1]], border = color_assign[label_data[sector.name,1]], ytop.column = 1, ybottom = 0, ...)  
        circos.lines(c(0, max(region)), c(0.5, 0.5), col = 'gray', lwd = 0.3)  
    } )

    category_legend <- Legend(
    labels = c(unique(label_data)),#各二级分类的名称
    type = 'points', pch = NA, background = c(outcolor), #各二级分类的颜色
    labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'),
    title = 'category')
    updown_legend <- Legend(
    labels = c('Up-regulated', 'Down-regulated'), 
    type = 'points', pch = NA, background = c(upregulated, downregulated), 
    labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'),
    title='differential expression')
    pvalue_legend <- Legend(
    col_fun = colorRamp2(round(seq(0, p_max, length.out = 6), 0), 
                        colorRampPalette(c(lengend_down, lengend_up))(6)),
    legend_height = unit(3, 'cm'), labels_gp = gpar(fontsize = 8), 
    title_gp = gpar(fontsize = 9), title_position = 'topleft', title = '-Log10(Pvalue)',
    direction = "horizontal")
    lgd_list_vertical <- packLegend(category_legend, updown_legend, pvalue_legend)
    grid.draw(lgd_list_vertical)
    # circos.clear()
}
