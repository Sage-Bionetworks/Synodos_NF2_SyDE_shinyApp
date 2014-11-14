
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)



shinyServer(function(input, output, session) {
  
  selected_pathways <- reactive({
    #selected grouped pathways
    selected_grps <- names(grouped_pathways) %in% input$grouped_pathways
    pathways <- unlist(grouped_pathways[selected_grps])
    #selected pathways based on entered genes
    pathways_based_on_selected_genes <- genes_to_pathways$pathway[genes_to_pathways$genes %in% selected_genes()]
    #now add any individual pathways
    pathways <- unique(c(pathways, input$pathways, pathways_based_on_selected_genes))
    pathways
  })
  
  
  genes_in_pathways <- reactive({
    #selected grouped pathways
    selected_grps <- names(grouped_pathways) %in% input$grouped_pathways
    pathways <- unlist(grouped_pathways[selected_grps])
    #now add any individual pathways
    pathways <- unique(c(pathways, input$pathways))
    #find the genes involved in all the selected pathways
    selected_genes_based_on_pathways <- genes_to_pathways[genes_to_pathways$pathway %in% pathways,][['genes']]
    selected_genes_based_on_pathways
  })  
  
  selected_genes <- reactive({
    input$searchByGenes_button
    selected_genes <- isolate(input$custom_gene_list)
    geneList <- unlist(strsplit(selected_genes,split=c('[\\s+,\\n+\\r+)]'),perl=T))
    #convert everything to upper case
    geneList <- toupper(geneList)
    geneList <- geneList[ !geneList == "" ] #remove the blank entries
    geneList
  })
  
  
  selected_samples <- reactive({
    selected_samples <- NULL
     if(! is.null(input$selected_geo_studies)){
       publicData_phenotype_flt <- filter(publicData_phenotype, studyId %in% input$selected_geo_studies)
       selected_samples <- unique(paste0(publicData_phenotype_flt$studyId, "_", publicData_phenotype_flt$sampleId))
     }
   selected_samples
  })

  
  filtered_pathway_enrichment_matrix <- reactive({
      validate(
        need(length(selected_pathways()) > 0, "Please select atleast one choice in option 1,2,3")
      )
      #filter based on pathways
      rows_to_keep <- rownames(pathway_enrichment_scores) %in% selected_pathways()
      filtered_matrix <- pathway_enrichment_scores[rows_to_keep,] 
      
      #filter based on selected studies
      if(! is.null(selected_samples())){
        rows_to_keep <- colnames(filtered_matrix) %in%  selected_samples()
        filtered_matrix <- filtered_matrix[,rows_to_keep]
      }
      filtered_matrix
  })
  
  
  filtered_expression_matrix <- reactive({
    #filter based on genes
    genes <- as.character(unique(c(selected_genes(), genes_in_pathways())))
    validate(
      need(length(genes) > 0, "Please select atleast one choice in option 1,2,3")
    )
    rows_to_keep <- rownames(publicData_expMat) %in% genes
    filtered_expression_matrix <- publicData_expMat[rows_to_keep,]  
        
    #filter based on selected studies
    if(! is.null(selected_samples())){
      rows_to_keep <- colnames(filtered_expression_matrix) %in%  selected_samples()
      filtered_expression_matrix <- filtered_expression_matrix[,rows_to_keep]
    }

    #filter rows with NA values
    filtered_expression_matrix <- na.omit(filtered_expression_matrix)
    
    #make sure the filtered data only has < 1000 genes
    validate(
      need(nrow(filtered_expression_matrix) < 1000, 
           paste0("Due to longer compute time Expression clustering will only work upto 1000 
                  genes \n currently entered ",nrow(filtered_expression_matrix), " genes")
      )
    )
    filtered_expression_matrix  
  })
  
  
  pubData_filteredAnnotation <- reactive({
    annotation <- publicData_phenotype[,c('studyId', 'sampleId')]
    annotation['studyType'] = 'meningioma'
    annotation$studyType[annotation$studyId %in% get_schwannoma_studies()] = 'schwannoma'
    
    #if only one geo study is selected selected 
    if( length(input$selected_geo_studies) == 1 ){
      annotation <- filter(publicData_phenotype, studyId %in% input$selected_geo_studies )
      annotation <- dcast(annotation, sampleId + studyId  ~ variable)
      cols_to_keep <- c('sampleId', 'studyId', input$pubData_selectedPhenotypes)
      annotation <- annotation[,cols_to_keep]
    }  
    #remove duplicated rows
    rows_to_keep <- ! duplicated(annotation)
    annotation <- annotation[rows_to_keep,]
    row.names(annotation) <- paste0(annotation$studyId, '_', annotation$sampleId)
    annotation$sampleId <- NULL
    annotation
  })
  
  
  
  get_schwannoma_studies <- reactive({
    schwannoma_studies
  })
  
  #heatmap for pubData pathway enrichment 
  output$pubData_pathway_heatmap <- renderPlot({
    fontsize_row=10
    if(nrow(filtered_pathway_enrichment_matrix()) > 100){
      fontsize_row = 0
    }
    scaled_mat = t(scale(t(filtered_pathway_enrichment_matrix())))
    withProgress(session, {
                            setProgress(message = "clustering & rendering heatmap, please wait", 
                                        detail = "This may take a few moments...")
                            expHeatMap(scaled_mat, pubData_filteredAnnotation(),
                                       scale=F,
                                       fontsize_col=0,
                                       fontsize_row=fontsize_row)
    })
  })
  
  #heatmap for EXPRESSION pubData 
  output$pubData_expression_heatmap <- renderPlot({  
    fontsize_row=10
    if(nrow(filtered_expression_matrix()) > 100){
      fontsize_row = 0
    }
    scaled_expmat = t(scale(t(filtered_expression_matrix())))
    withProgress(session, {
      setProgress(message = "clustering & rendering heatmap, please wait", 
                  detail = "This may take a few moments...")
      expHeatMap(scaled_expmat, 
                 pubData_filteredAnnotation(),
                 scale = F,
                 fontsize_col=0,
                 fontsize_row=fontsize_row)
    })
  })
 
  #update the available phenotypes for a single study
  observe({
    phenotypes <- filter(publicData_phenotype, studyId %in% input$selected_geo_studies)
    choices = unique(as.character(phenotypes[['variable']]))
    updateSelectInput(session = session,
                      inputId = "pubData_selectedPhenotypes",
                      choices = choices)
  })

  ###############################
  ## Kinome Screen Section
  ###############################
  output$kinome_barPlot <- renderChart({
    p <- nPlot(log2_ratio ~ Uniprot , data=get_ordered_kinomeData(), group="condition", 
               type="multiBarChart")
    p$params$width <- 800
    p$params$height <- 500
    p$addParams(dom = 'kinome_barPlot')
    return(p)
  })
  
  get_cleaned_kinomeData <- reactive({
    filtered_data <- kinomeData[!apply(kinomeData, 1, function(x) {sum(is.na(x)) > 0}),]
    filtered_data
  })
  
  get_missing_kinomeData <- reactive({
    missing_data <- finalData[apply(kinomeData, 1, function(x) {sum(is.na(x)) > 0}),]
    missing_data
  })
  
  get_filtered_kinomeData <- reactive({
     cleaned_data <- get_cleaned_kinomeData()
     #filter kinase family
     if( sum(input$kinome_selected_kinaseFamily %in% c('ALL')) != 1 ){
       cleaned_data <- filter(cleaned_data, Family %in% input$kinome_selected_kinaseFamily )        
     }
     #filter by samples
     cleaned_data <- filter(cleaned_data, condition %in% input$kinome_selected_samples )        
     cleaned_data
   })

  get_ordered_kinomeData <- reactive({
    kinomeData <- get_filtered_kinomeData()
    #get the mean of ratios per gene across all the samples
    new_order <- names(sort(tapply(kinomeData$log2_ratio, kinomeData$Gene, mean), decreasing=T))
    
    new_order <- data.frame('order'= match(kinomeData$Gene,new_order),'index' = 1:nrow(kinomeData)) %>%
                   arrange(order) %>%
                   select(index)
    kinomeData[new_order$index,]
  })

  ##########################
  # Drug Screening Section
  ##########################
  
  get_drug_flt_normViab <- reactive({
    flt_Drug_normViab <- drug_normViab
    if(! is.null(input$selected_drugs)){
      flt_Drug_normViab <- filter(flt_Drug_normViab, drug %in% input$selected_drugs)  
    }
    if(! is.null(input$drugs_to_remove)){
      flt_Drug_normViab <- filter(flt_Drug_normViab, ! drug %in% input$drugs_to_remove)  
    }
    if(! is.null(input$selected_cellLines)){
      flt_Drug_normViab <- filter(flt_Drug_normViab, cellLine %in% input$selected_cellLines)  
    }
    flt_Drug_normViab
  })
  
  get_drug_flt_ICVals <- reactive({
    flt_drug_ICVals <- drug_ICVals
    if(! is.null(input$selected_drugs)){
      flt_drug_ICVals <- filter(flt_drug_ICVals, drug %in% input$selected_drugs)  
    }
    if(! is.null(input$drugs_to_remove)){
      flt_drug_ICVals <- filter(flt_drug_ICVals, ! drug %in% input$drugs_to_remove)  
    }
    if(! is.null(input$selected_cellLines)){
      flt_drug_ICVals <- filter(flt_drug_ICVals, cellLine %in% input$selected_cellLines)  
    }
    flt_drug_ICVals
  })
  
  output$drugScreen_ICx_plot <- renderPlot({
    flt_drug_ICVals <- get_drug_flt_ICVals()
    ICx <- eval(paste0('IC', input$selected_IC_value))
    flt_drug_ICVals[ICx] <- log10(as.numeric(flt_drug_ICVals[,ICx]))  
    drug_levels <- flt_drug_ICVals %>%
                      group_by(drug) %>%
                      summarise(mean=mean(IC50, na.rm=T)) %>%
                      arrange(desc(mean)) %>% select(drug)
    drug_levels <- drug_levels$drug
    flt_drug_ICVals$drug <- factor(flt_drug_ICVals$drug,levels=drug_levels)
    facet_by <- paste(input$facet_by, collapse = ' + ' )
    facet_by <- formula(paste(facet_by, ' ~ .'))
    p <- ggplot(data=flt_drug_ICVals, aes_string(x="drug", y=ICx, group="cellLine")) + geom_line(aes(color=cellLine)) 
    p <- p +  geom_point(color='grey50') + theme_bw()
    if( length(input$facet_by) > 0){
       p <- p  + facet_grid(facet_by)    
    }
    p + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab('Drug') + ylab(paste0(ICx, ' (log10 molar conc)'))
  })

  
  
  
  output$drugResponse_plots <- renderPlot({
    validate(need(length(input$selected_drugs) != 0, paste0(" Please select drug/s (max upto 4)")))  
    validate(need(length(input$selected_drugs) < 5, paste0(" Please select < 5 drugs")))  
    flt_drug_normViab <- get_drug_flt_normViab()
  
    
    doseResp <- ddply(.data=flt_drug_normViab, .variables = c('drug', 'cellLine', 'experiment', 'group'), 
                      .fun = tmp_iterator, .parallel = T)
    facet_by <- paste(input$facet_by, collapse = ' + ' )
    facet_by <- formula(paste(facet_by, ' ~ drug'))
    p <- ggplot(data=doseResp, aes(x=fittedX, y=fittedY, group=cellLine))
    p <- p + geom_line(aes(color=cellLine)) + facet_grid( facet_by) + theme_bw()
    p <- p + geom_hline(aes(yintercept=0.5), color='grey50', linetype='dashed')
    p <- p + xlab('molar conc (log10)') + ylab('cell viability %')
    p
  })

})



# ########################
# ## TEST SECTION
# ########################
# 
#   ##test rchart
#   output$chart1 <- renderChart({
#     hmap <- rPlot(variable ~ School, color = 'rescale', data = findatamelt, type = 'tile')
#     hmap$addParams(height = 400, width=1000, dom="chart1")
#     hmap$guides(reduceXTicks = FALSE)
#     hmap$guides("{color: {scale: {type: gradient, lower: white, upper: red}}}")
#     hmap$guides(y = list(numticks = length(unique(findatamelt$value))))
#     hmap$guides(x = list(numticks = 5))
#     return(hmap)
#   })
# 
#   output$myChart <- renderChart({
#     names(iris) = gsub("\\.", "", names(iris))
#     p1 <- rPlot(input$x, input$y, data = iris, color = "Species", 
#                 facet = "Species", type = 'point')
#     p1$addParams(dom = 'myChart')
#     return(p1)
#   })

