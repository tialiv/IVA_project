# Volcano plot for multiple group

d1 <- final_all
top <- list()
for (i in 1:12) {
  d1[[i]][is.na(d1[[i]]$padj),] <- 1
  d1[[i]]$label <- ifelse(d1[[i]]$padj < 0.1, "FDR<0.1", "FDR>0.1")
  top[[names(d1[i])]] <- top_n(d1[[i]], 20, abs(d1[[i]]$log2FoldChange))
}

top_gene <- data.frame()
for (i in 1:12) {
  temp <- top[[i]] %>% mutate(group = rep(names(top[i]), nrow(top[[i]])))
  top_gene <- rbind(top_gene, temp)
}

all <- data.frame()
for (i in 1:12) {
  temp <- d1[[i]] %>% mutate(group = rep(names(d1[i]), nrow(d1[[i]])))
  all <- rbind(all, temp)
}

all$size <- case_when(!(all$geneid %in% top_gene$geneid) ~ 1, 
                      all$geneid %in% top_gene$geneid ~ 2)

dt <- filter(all, size == 1)
theme <- theme(text=element_text(size=15),
               axis.text.x=element_text(color="black", angle = 90, vjust = 0.3, hjust = 0.3),
               axis.text.y=element_text(color="black"),
               axis.ticks.x.bottom=element_blank(),
               #axis.ticks.length.x=unit(-.25, "cm"),
               strip.background = element_rect(colour=NA, fill=NA),
               plot.background = element_blank(),
               #panel.grid.major = element_blank(),
               #panel.grid.minor = element_blank(),
               panel.spacing.x=unit(.5, "lines"),
               panel.border = element_rect(color = "black", fill = NA, size = 1),
               axis.line.x = element_line(size=0.4),
               axis.line.y = element_line(size=0.4),
               legend.key.size = unit(2, 'lines'),
               panel.background = element_blank()
)

dfcol <- data.frame(x=c(1:12), y = 0, label = unique(all$group))
mycol <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", 
  "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f")
ggplot() + theme +
  geom_jitter(data = dt, mapping = aes(x=group, y=log2FoldChange, color = label), size = 0.85, width = 0.4,
              alpha = 0.5) +
  geom_jitter(data = top_gene, mapping = aes(x=group, y=log2FoldChange, color = label)) + ylim(c(-15,15)) +
  #geom_tile(data = dfcol, aes(x=x, y=y), height = 5, color = "black", fill = mycol, alpha = 0.6, show.legend = F) +
  geom_text_repel(data = top_gene, aes(x= group, y = log2FoldChange, label = external_gene_name),
                  force = 1.2, arrow = arrow(length = unit(0.008, "npc"),
                                             type = "open", ends = "last")) +
  scale_color_manual(name = NULL, values = c("red", "black")) +
  labs(x="Group", y = "Log2 Fold Change") #+ 
  #geom_text(data = dfcol, aes(x, y, label = label), color = "white", size = 1)
  
  
## Interactive app for DESeq2 results
library(DT)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(shinydashboard)
library(circlize)
library(GetoptLong)
library(shiny)

dds <- dataset[["dds_DES"]]
res <- comparison01[[1]]

res <- comparison01
res <- get_id(res)

env = new.env()

make_heatmap = function(fdr = 0.01, base_mean = 0, log2fc = 1) {
  l = res$padj <= fdr & res$baseMean >= base_mean & 
    abs(res$log2FoldChange) >= log2fc; l[is.na(l)] = FALSE
    
    if(sum(l) == 0) return(NULL)
    
  #  re <- res[l,]
    m = counts(dds, normalized = TRUE) %>% 
      data.frame() %>%
      tibble::rownames_to_column(var="gene") %>%
      as_tibble() %>%
      left_join(output, by=c("gene" = "ensembl_gene_id")) 
    id = match(res$geneid, m$gene)
    m = m[id, ]
    m = m[l,]
    
    env$row_index = which(l)
    
    ht = Heatmap(t(scale(t(m[,2:(ncol(m)-4)]))), name = "z-score",
                 top_annotation = HeatmapAnnotation(
                   culture = colData(dds)$culture, 
                   col = list(culture = c("Fresh" = "#ffff99", "Cultured" = "#beaed4")),
                   sizeFactor = anno_points(colData(dds)$sizeFactor)
                 ),
                 show_row_names = FALSE, show_column_names = FALSE, row_km = 2,
                 column_title = paste0(sum(l), " significant genes with FDR < ", fdr),
                 show_row_dend = FALSE) + 
      Heatmap(log10(res$baseMean[l]+1), show_row_names = FALSE, width = unit(5, "mm"),
              name = "log10(baseMean+1)", show_column_names = FALSE) +
      Heatmap(res$log2FoldChange[l], show_row_names = FALSE, width = unit(5, "mm"),
              name = "log2FoldChange", show_column_names = FALSE,
              col = colorRamp2(c(-2, 0, 2), c("#377eb8", "white", "#fc8d62")))
    ht = draw(ht, merge_legend = TRUE)
    ht
}

# make the volcano plot with some genes highlited
make_volcano = function(res, highlight = NULL) {
  col = rep("#00000020", nrow(res))
  cex = rep(0.5, nrow(res))
  names(col) = res$external_gene_name
  names(cex) = res$external_gene_name
  if(!is.null(highlight)) {
    col[highlight] = "red"
    cex[highlight] = 1
  }
  x = res$log2FoldChange
  y = -log10(res$padj)
  col[col == "red" & x < 0] = "darkgreen"
  par(mar = c(4, 4, 1, 1))
  
  suppressWarnings(
    plot(x, y, col = col, 
         pch = 16, 
         cex = cex,
         xlab = "log2 fold change", ylab = "-log10(FDR)")
  )
}

body = dashboardBody(
  fluidRow(
    column(width = 4,
           box(title = "Differential heatmap", width = NULL, solidHeader = TRUE, status = "primary",
               originalHeatmapOutput("ht", height = 800, containment = TRUE)
           )
    ),
    column(width = 4,
           id = "column2",
           box(title = "Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary",
               subHeatmapOutput("ht", title = NULL, containment = TRUE)
           ),
           box(title = "Output", width = NULL, solidHeader = TRUE, status = "primary",
               HeatmapInfoOutput("ht", title = NULL)
           ),
           box(title = "Note", width = NULL, solidHeader = TRUE, status = "primary",
               htmlOutput("note")
           ),
    ),
    column(width = 3,
           box(title = "Volcanno plot", width = NULL, solidHeader = TRUE, status = "primary",
               plotOutput("volcanno_plot")
           ),
           box(title = "Result table of the selected genes", width = NULL, solidHeader = TRUE, status = "primary",
               DTOutput("res_table")
           )
    ),
    tags$style("
            .content-wrapper, .right-side {
                overflow-x: auto;
            }
            .content {
                min-width:1500px;
            }
        ")
  )
)

brush_action = function(df, input, output, session) {
  
  row_index = unique(unlist(df$row_index))
  selected = env$row_index[row_index]
  
  output[["volcanno_plot"]] = renderPlot({
    make_volcano(res, selected)
  })
  
  output[["res_table"]] = renderDT(
    formatRound(datatable(res[selected, c("baseMean", "log2FoldChange", "padj", "external_gene_name")], 
                          rownames = F), columns = 1:3, digits = 3)
  )
  
  output[["note"]] = renderUI({
    if(!is.null(df)) {
      HTML(qq("<p>Row indices captured in <b>Output</b> only correspond to the matrix of the differential genes. To get the row indices in the original matrix, you need to perform:</p>
<pre>
l = res$padj <= @{input$fdr} & 
    res$baseMean >= @{input$base_mean} & 
    abs(res$log2FoldChange) >= @{input$log2fc}
l[is.na(l)] = FALSE
which(l)[row_index]
</pre>
<p>where <code>res</code> is the complete data frame from DESeq2 analysis and <code>row_index</code> is the <code>row_index</code> column captured from the code in <b>Output</b>.</p>"))
    }
  })
}

ui = dashboardPage(
  dashboardHeader(title = "DESeq2 results"),
  dashboardSidebar(
    selectInput("fdr", label = "Cutoff for FDRs:", c("0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05, "0.1" = 0.1)),
    numericInput("base_mean", label = "Minimal base mean:", value = 0),
    numericInput("log2fc", label = "Minimal abs(log2 fold change):", value = 0),
    actionButton("filter", label = "Generate heatmap")
  ),
  body
)

server = function(input, output, session) {
  observeEvent(input$filter, {
    ht = make_heatmap(fdr = as.numeric(input$fdr), base_mean = input$base_mean, log2fc = input$log2fc)
    if(!is.null(ht)) {
      makeInteractiveComplexHeatmap(input, output, session, ht, "ht",
                                    brush_action = brush_action)
    } else {
      # The ID for the heatmap plot is encoded as @{heatmap_id}_heatmap, thus, it is ht_heatmap here.
      output$ht_heatmap = renderPlot({
        grid.newpage()
        grid.text("No row exists after filtering.")
      })
    }
  }, ignoreNULL = FALSE)
}

shinyApp(ui, server)

interactivate(dds)

body = dashboardBody(
  fluidPage(
    titlePanel("DESeq2 results visualization"),
      mainPanel(tabsetPanel(
        tabPanel("DESeq2 results",
  fluidRow(
    column(width = 3,
           box(title = "Differential heatmap", width = NULL, solidHeader = TRUE, status = "primary",
               originalHeatmapOutput("ht", height = 800, containment = TRUE)
           )
    ),
    column(width = 3,
           id = "column2",
           box(title = "Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary",
               subHeatmapOutput("ht", title = NULL, containment = TRUE)
           ),
           box(title = "Output", width = NULL, solidHeader = TRUE, status = "primary",
               HeatmapInfoOutput("ht", title = NULL)
           ),
           box(title = "Note", width = NULL, solidHeader = TRUE, status = "primary",
               htmlOutput("note")
           ),
    ),
    column(width = 3,
           box(title = "Volcanno plot", width = NULL, solidHeader = TRUE, status = "primary",
               plotOutput("volcanno_plot")
           ),
           box(title = "Result table of the selected genes", width = NULL, solidHeader = TRUE, status = "primary",
               DTOutput("res_table")
           )
    ))),
    tabPanel("Functional analysis results",
             fluidRow(
               column(width = 4,
                      box(title = "Result table of the selected genes", width = NULL, solidHeader = TRUE, status = "primary",
                          DTOutput("res_table")
                      )
               ),
               column(width = 4,
                      id = "column2",
                      box(title = "GO enrichment", width = NULL, solidHeader = TRUE, status = "primary",
                          subHeatmapOutput("ht", title = NULL, containment = TRUE)
                      )
               )),
    tags$style("
            .content-wrapper, .right-side {
                overflow-x: auto;
            }
            .content {
                min-width:1500px;
            }
        ")
      )
    )
   )
  )
)


