kinomeViewerModuleUI <- function(id){
  ns <- NS(id)
  tagList(
  myHeader <- dashboardHeader(disable=TRUE),
  mySidebar <- dashboardSidebar(disable=TRUE),
  
  myBody <- dashboardBody(
    tabBox(width = 12,
             tabPanel("Baseline Differential Expression",
                      column(width = 3,style='padding-right:0px;',
                             box(width=NULL, status='primary', collapsible=FALSE, 
                                 collapsed=FALSE, solidHeader=TRUE,
                                 title = tagList(shiny::icon("check", lib="glyphicon"),
                                                 "Select genes"),
                                 uiOutput(ns("genesbase")))),
                      column(width = 8,style='padding-left:30px;',
                             plotOutput(ns("waterfall.hs"), height = 300),
                             plotOutput(ns("waterfall.syn"), height = 300))),
      tabPanel("Treatment Heatmap",
    tags$head(tags$style(HTML('
                                .col-sm-3, {
                                  padding:0 !important;
                                }
                              '))),
    fluidRow(
      column(width = 3,style='padding-right:0px;',
             
             # Choose sample labels
             box(width=NULL, status='primary', collapsible=TRUE, 
                 collapsed=FALSE, solidHeader=TRUE,
                 title = tagList(shiny::icon("th-list", lib="glyphicon"),
                                 "Label samples"),               
                 uiOutput(ns("anno"))
                
             ),
             
             #kinase input
             box(width=NULL, status='primary', collapsible=TRUE, 
                 collapsed=FALSE, solidHeader=TRUE,
                 title = tagList(shiny::icon("th-list", lib="glyphicon"),
                                 "Kinase Input"),               
                 uiOutput(ns("genes"))
             ),
             
             #Clustering box
             box(width = NULL, status = "warning", solidHeader=TRUE, 
                 collapsible=TRUE, collapsed=TRUE,
                 title = tagList(shiny::icon("wrench", lib="glyphicon"), "Change cluster options"),
                 #distance metric
                 selectInput(ns("clustering_distance"), "Distance Calculation",
                             choices=c("correlation", "euclidean", "maximum", 
                                       "manhattan", "canberra", "binary", "minkowski"),
                             selectize=T, multiple=F, selected="euclidean"),
                 # set the clustering method
                 selectInput(ns("clustering_method"), "Clustering Method",
                             choices=c("ward", "single", "complete", "average", 
                                       "mcquitty", "median", "centroid"),
                             selectize=T, multiple=F, selected="average"),
                 checkboxInput(ns('cluster_cols'), 'Cluster the columns', value = TRUE),
                 checkboxInput(ns('cluster_rows'), 'Cluster the rows', value = TRUE)
             )
     ),
     column(width = 8,style='padding-left:30px;',
            plotOutput(ns("heatmap"), height = 650)
                          )
             )
    )
  )
  ))
  dashboardPage(header=myHeader, sidebar=mySidebar, body=myBody,
                skin = "blue")
}

kinomeViewerModule <- function(input,output,session,data,tag){
  dataset <- reactive({
    ds <- data
    flog.debug(sprintf("filtered ds dims: %s", dim(ds)), name="server")
    #rows_to_keep <- apply(exprs(ds), 1, var) > 0.1
    rows_to_keep <- order(apply(exprs(ds),1,var),decreasing=T)
    ds_filtered <- ds[rows_to_keep, ]
    
    ds_filtered
  })
  
  metadata <- reactive({
    m_eset <- dataset()
    metaData <- pData(m_eset)
    metaData
  })

  ns <- NS(tag)
  output$anno <- renderUI({
    metaData <- metadata()
    tagList(
      selectInput(ns('annotation_labels'),'Annotate Samples by:',
                  choices=colnames(metaData), selectize=T, multiple=T, selected=colnames(metaData)[1])
    )
  })
  
  
  output$genes <- renderUI({
    m_eset <- dataset()
    rows_to_keep <- order(apply(exprs(m_eset),1,var),decreasing=T)
    m_top500 <- m_eset[rows_to_keep,]
    geneList <- rownames(m_top500)
    geneList <- na.omit(geneList)
    tagList(
    tags$textarea(paste0(c(geneList), collapse="\n"), rows=5, id=ns("selected_genes"), style="width: 100%"),
    actionButton(ns("refreshGene"), "Refresh")
    )
  })
  
  output$genesbase <- renderUI({
    m_eset <- dataset()
    geneList <- rownames(m_eset)
    geneList <- na.omit(geneList)
    tagList(
      tags$textarea(paste0(c(geneList), collapse="\n"), rows=5, id=ns("selected_genes2"), style="width: 100%"),
      actionButton(ns("refreshGene2"), "Refresh")
    )
  })

  
  user_submitted_selections <- reactive({
    ds <- dataset()
      input$refreshGene
        geneList <- isolate(input$selected_genes)
        geneList <- clean_list(geneList)
        geneList<- intersect(geneList, rownames(fData(ds)))
    geneList
  })

  filtered_dataset <- reactive({
    ds <- dataset()
    selected_genes <- user_submitted_selections()
    validate(need(length(selected_genes) > 4, "Please select at least 5 genes." ))
    validate(need(length(selected_genes) < 2001, "Number of total selected genes is at most 2000."))

    ds <- ds[selected_genes,]

    ds
  })
  
  heatmap_cache <- reactiveValues()
  
  anno_labels <- reactive({
    validate(need(length(input$annotation_labels) <= 2, "Please select at most 2 labels."))
    input$annotation_labels
  })
  
  #return the heatmap plot
  output$heatmap <- renderPlot({  
    flog.debug("Making heatmap", name='server')
    
    cluster_rows <- input$cluster_rows
    cluster_cols <- input$cluster_cols
    
    m_eset <- filtered_dataset()
    m <- exprs(m_eset)
    
    keep <- rowSums(is.na(m)) < 3
    m <- m[keep, ] 
    
    m <- data.matrix(m)
    
    validate( need( ncol(m) != 0, "Filtered matrix contains 0 samples.") )
    validate( need( nrow(m) != 0, "Filtered matrix contains 0 genes.") )
    
    metadata <- metadata()
    anno <- anno_labels()
    annotation <- get_heatmapAnnotation(anno, metadata)
    
    fontsize_row <- ifelse(nrow(m) > 100, 0, 8)
    fontsize_col <- ifelse(ncol(m) > 50, 0, 8)    
    
    heatmap.color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
    
    
    heatmap_cache$heatmap <- sbHeatMap(m,annotation,
                                        clustering_distance_rows = input$clustering_distance,
                                        clustering_distance_cols = input$clustering_distance,
                                        fontsize_col=fontsize_col, 
                                        fontsize_row=fontsize_row,
                                        scale=F,
                                        color=heatmap.color,
                                        #breaks=heatmap.breaks,
                                        clustering_method = input$clustering_method,
                                        explicit_rownames = fData(m_eset)$explicit_rownames,
                                        cluster_rows=cluster_rows, cluster_cols=cluster_cols,
                                        drawColD=FALSE)
    
  })
  
  output$waterfall.hs <- renderPlot({
    
    user_submitted_selections <- reactive({
      ds <- dataset()
      input$refreshGene2
      geneList <- isolate(input$selected_genes2)
      geneList <- clean_list(geneList)
      geneList<- intersect(geneList, rownames(fData(ds)))
      geneList
    })
    
    geneList <- user_submitted_selections()
    hs.filt <- HS01.HS11.base %>% filter(Gene %in% c(geneList, "ULK4", "EPHA3", "CAMK4", "STK17B", "CHEK1", "AURKA"))
    p<-ggplot(data = hs.filt, aes(x=Gene, 
                                      y=mean.log2ratio, 
                                      fill = mean.log2ratio, group = comp)) +
      geom_errorbar(aes(x=Gene, ymin=mean.log2ratio-sem, ymax = mean.log2ratio+sem), stat = "identity", position = "dodge") +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_viridis(option="plasma") +
      scale_color_viridis(option="plasma") +
      coord_flip() +
      labs(y = "log2ratio (mean)", x = "Gene", title = "HS01 vs HS11, baseline") + 
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
    
    p
  })
    
    output$waterfall.syn <- renderPlot({
    
      user_submitted_selections <- reactive({
        ds <- dataset()
        input$refreshGene2
        geneList <- isolate(input$selected_genes2)
        geneList <- clean_list(geneList)
        geneList<- intersect(geneList, rownames(fData(ds)))
        geneList
      })
      
      geneList <- user_submitted_selections()
      syn.filt <- Syn5.Syn1.base %>% filter(Gene %in% c(geneList, "EPHB1", "EPHA4", "KIT", "AURKA", "MAP3K10", "STYK1"))

    p<-ggplot(data = syn.filt, aes(x=Gene, 
                                      y=mean.log2ratio, 
                                      fill = mean.log2ratio, group = comp)) +
      geom_errorbar(aes(x=Gene, ymin=mean.log2ratio-sem, ymax = mean.log2ratio+sem),stat = "identity", position = "dodge") +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_viridis(option="plasma") +
      scale_color_viridis(option="plasma") +
      coord_flip() +
      labs(y = "log2ratio (mean)", x = "Gene", title = "Syn5 vs Syn1, baseline") + 
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
    
    p
  })

}


grid_arrange_shared_legend <- function(plots, ncol = length(plots), nrow = 1, position = c("bottom", "right")) {  
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}

clean_list <- function(x) {
  # Split by space, comma or new lines
  x <- as.character(x)
  x <- unlist(strsplit(x, split=c('[\\s+,\\n+\\r+)]'),perl=T))
  
  # remove the blank entries
  x <- x[!(x == "")]
  
  x
}
