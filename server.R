
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyBS)


shinyServer(function(input, output, session) {
  

  ###############################
  ## Kinome Screen Section
  ###############################
  
  #kinome_selected_points <- linked_brush(keys=NULL, "red")
  
  #two part filtering is needed to make sure histograms are plotted using data which is 
  # only filtered by samples and kinases and NOT by user selected thresholds
  #filer 1 : 
  get_filtered_kinomeData_by_Samples_N_Kinases <- reactive ({
    filtered_data <- kinomeData
    #1. filter by samples
    validate( need( length(input$kinome_selected_samples) != 0, "Please select atleast one sample from option 1"))
    filtered_data <- filter(filtered_data, condition %in% input$kinome_selected_samples )        
    
    #2. filter kinase family  and selected genes
    num_selected_entities <- length(input$kinome_selected_kinaseFamily) + length(input$kinome_selected_genes)
    validate( need( num_selected_entities != 0, "Please select atleast one kinase family and/or gene/s"))
    if( sum(input$kinome_selected_kinaseFamily %in% c('ALL')) != 1 ){
      filtered_data <- filter(filtered_data, Family %in% input$kinome_selected_kinaseFamily | Gene %in% input$kinome_selected_genes)        
    }
    
    #dynamically update sliders to reflect max values based on current selection
    updateSliderInput(session, "kinome_var_threshold", max = max(filtered_data$variability, na.rm=T)+1,
                      value = c(0,max(filtered_data$variability, na.rm=T)+1 ))
    updateSliderInput(session, "kinome_PeptideCount_threshold", max = max(filtered_data$count, na.rm=T),
                      value = c(0,max(filtered_data$count, na.rm=T) ))
    return(filtered_data)
  })

  #filer 2: now filtering based on user selected thresholds
  get_filtered_kinomeData_by_userSelected_thresholds <- reactive({
      filtered_data <-  get_filtered_kinomeData_by_Samples_N_Kinases()
      
      #A. filter based on % variability
      min_var <- input$kinome_var_threshold[1]
      max_var <- input$kinome_var_threshold[2]
      # keeping those kinases which have variability = NA by 'is.na(variability)'
      # we dont want to remove these 
      # the variability is NA as the peptide count  == 1 
      filtered_data <- filter(filtered_data, (variability >= min_var & variability <= max_var) | is.na(variability) )
          
      #B.  filter based on peptide counts
      min_count_threshold <- input$kinome_PeptideCount_threshold[1]
      max_count_threshold <- input$kinome_PeptideCount_threshold[2]
      filtered_data <- filter(filtered_data, count >= min_count_threshold & count <= max_count_threshold)
      
      #C.filter based on ratio
      #only include those kinases which fall between the given ratios
      min_ratio <- input$kinome_ratio_includeRegion[1]
      max_ratio <- input$kinome_ratio_includeRegion[2]
      filtered_data <- filter(filtered_data, abs(log2ratio) >= log2(min_ratio) & abs(log2ratio) <= log2(max_ratio)) 
      
      #create a id columns for tooltip
      filtered_data['id'] <- get_tooltip_vals(filtered_data)
      
      #kinome_selected_points$set_keys(seq_len(nrow(filtered_data)))
      return(filtered_data)
  })
  
  get_tooltip_vals <- function(df){
    selected_cols <- c('Description', 'Gene', 'Uniprot', 'Family',
                       'ratio', 'variability', 'count', 'uniq_peptides')
    df <- df[,selected_cols]
    apply(df,1, function(x) paste0(names(x), ": ", format(x), collapse = "<br />"))
  }
  
  
  get_ordered_kinomeData <- reactive({
    kinomeData <- get_filtered_kinomeData_by_userSelected_thresholds()
    #get the mean of ratios per gene across all the samples
    new_order <- names(sort(tapply(kinomeData$log2ratio, kinomeData$Gene, mean), decreasing=T))
    new_order <- data.frame('order'= match(kinomeData$Gene,new_order),'index' = 1:nrow(kinomeData)) %>%
                  arrange(order) %>%
                  select(index)
    kinomeData[new_order$index,]
  })
    
  #custom height deciding function for kinome barplot 1
  get_plotHeight <- reactive({
    x <- get_filtered_kinomeData_by_userSelected_thresholds()
    if( nrow(x) < 40 )  return(400)
    else nrow(x) * 10  
  })
  
  output$kinome_barPlot <- renderPlot({
    withProgress(message = "computing", 
                 detail = "This may take a few moments...", {
                  x <-  get_filtered_kinomeData_by_userSelected_thresholds ()
                  new_levels <- unique(x[order(x$log2ratio),'Gene'])
                  x$Gene <- factor(x$Gene, levels=new_levels)
                  p <- ggplot(data=x, aes(y=log2ratio, x=Gene, fill=condition)) 
                  p <- p + geom_bar(stat="identity", position="dodge", width=0.5) 
                  p <- p + geom_errorbar(aes(ymax=log2ratio_max,ymin=log2ratio_min), width=0.5, position="dodge", colour="grey60") 
                  p <- p + coord_flip() + theme_bw() + theme(legend.title=element_blank(), legend.position="top")
                  p <- p + ylab('ratio(log2)')
                 })
  return(p)
  }, height=function(){get_plotHeight()})

  
  #ggvis interactive plot  
  get_filtered_kinomeData_by_userSelected_thresholds %>%
    ggvis(x = ~count, y = ~log2ratio, fill = ~condition, key := ~id, size.hover := 200) %>%
    layer_points() %>%
    add_legend(c("fill")) %>%
    add_tooltip(function(df) df$id, "hover") %>%
    set_options(width=550, height=350) %>%     
    bind_shiny("kinomeData_scatterPlot") 


  #3. Histogram of variation
  output$kinome_var_histogram <- renderPlot({
    x <- get_filtered_kinomeData_by_Samples_N_Kinases()
    p <- ggplot(data=x, aes(x=variability,fill=condition)) + geom_histogram(binwidth=5)
    p <- p + theme_bw() + xlab('percent variability') + ylab('') 
    p <- p + theme(legend.position="none")
    p <- p + geom_vline(xintercept = input$kinome_var_threshold, linetype='dashed', color="grey30")
    return(p)
  })

  #4. Histogram of peptide count variation
  output$kinome_uniqPeptides_histogram <- renderPlot({
   x <- get_filtered_kinomeData_by_Samples_N_Kinases()
    p <- ggplot(data=x, aes(x=count, fill=condition)) + geom_histogram(binwidth=1)
    p <- p + theme_bw() + xlab('#peptides counts') + ylab('')
    p <- p + theme(legend.position="none")
    p <- p + geom_vline(xintercept = input$kinome_PeptideCount_threshold, linetype='dashed', color="grey30")
    return(p)
  })

  #5. Histogram of ratio variation
  output$kinome_ratio_histogram <- renderPlot({
    x <- get_filtered_kinomeData_by_Samples_N_Kinases()
    p <- ggplot(data=x, aes(x=log2ratio, fill=condition)) + geom_histogram(binwidth=.05)
    p <- p + theme_bw() + xlab('ratio(log2)') + ylab('')
    p <- p + theme(legend.position="none")
    return(p)
  })


###################################################  
# Drug Screening Section
###################################################
  get_selected_cellLines <- reactive({
    MGH_cellLines <- if('ALL' %in% input$MGH_cellLines) unique(MGH_normViab$cellLine) else input$MGH_cellLines
    UCF_cellLines <- if('ALL' %in% input$UCF_cellLines) unique(UCF_normViab$cellLine) else input$UCF_cellLines
    cellLines <- c(MGH_cellLines,UCF_cellLines)
    validate(need(length(cellLines) != 0, "Atleast one cellLine needs to be selected" ) )
    cellLines
  })
  
  get_selected_drugs <- reactive({
    drugs <- if('ALL' %in% input$selected_drugs) unique(drug_normViab$drug) else input$selected_drugs
    validate(need(length(drugs) != 0, "Atleast one drug needs to be selected"))
    drugs
  })

  get_drug_flt_normViab <- reactive({
    flt_Drug_normViab <- filter(drug_normViab, drug %in% get_selected_drugs())  
    flt_Drug_normViab <- filter(flt_Drug_normViab, cellLine %in% get_selected_cellLines())  
    flt_Drug_normViab
  })
  
  get_drug_flt_ICVals <- reactive({
    flt_drug_ICVals <- filter(drug_ICVals, drug %in% get_selected_drugs())  
    flt_drug_ICVals <- filter(flt_drug_ICVals, cellLine %in% get_selected_cellLines())  
    flt_drug_ICVals
  })
  
  output$drugScreen_ICx_plot <- renderPlot({
    flt_drug_ICVals <- get_drug_flt_ICVals()
    ICx <- eval(paste0('IC', input$selected_IC_value))
    #remove NA
    flt_drug_ICVals <- flt_drug_ICVals[! is.na(flt_drug_ICVals[ICx]), ]
    #convert to log10
    flt_drug_ICVals[ICx] <- log10(as.numeric(flt_drug_ICVals[,ICx])) 
    #keep rows where log10 ICx <= 0
    flt_drug_ICVals <- flt_drug_ICVals[flt_drug_ICVals[ICx] <= 0,]
    
  
    drug_levels <- flt_drug_ICVals %>%
                      group_by(drug) %>%
                      summarise(med=median(IC50, na.rm=T)) %>%
                      arrange(desc(med)) %>% select(drug)
    drug_levels <- drug_levels$drug
    flt_drug_ICVals$drug <- factor(flt_drug_ICVals$drug,levels=drug_levels)
    facet_by <- paste(input$facet_by, collapse = ' + ' )
    facet_by <- formula(paste(facet_by, ' ~ .'))
    p <- ggplot(data=flt_drug_ICVals, aes_string(x="drug", y=ICx, group="cellLine")) 
    p <- p + geom_point(aes(color=cellLine), size=3) + theme_bw()
    if( length(input$facet_by) > 0){
       p <- p  + facet_grid(facet_by)    
    }
    p + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab('Drug') + ylab(paste0(ICx, ' (log10 molar conc)'))
  })
  

  output$drug_efficacy <- renderPlot({
    flt_drug_ICVals <- get_drug_flt_ICVals()
    
    drug_levels <- flt_drug_ICVals %>%
      group_by(drug) %>%
      summarise(med=median(maxEfficacy, na.rm=T)) %>%
      arrange(desc(med)) %>% select(drug)
    drug_levels <- drug_levels$drug
    flt_drug_ICVals$drug <- factor(flt_drug_ICVals$drug,levels=drug_levels)
    facet_by <- paste(input$facet_by, collapse = ' + ' )
    facet_by <- formula(paste(facet_by, ' ~ .'))
    p <- ggplot(data=flt_drug_ICVals, aes(x=drug, y=maxEfficacy*100, group=cellLine)) 
    p <- p + geom_point(aes(color=cellLine), size=3) + theme_bw()
    if( length(input$facet_by) > 0){
      p <- p  + facet_grid(facet_by)    
    }
    p <- p + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab('Drug') + ylab('% Efficacy')
    p
  })
  
  
  output$global_drugViab_heatMap <- renderPlot({
    
    validate(need(length(get_selected_cellLines()) != 0, paste0(" Please select cellLine/s")))  
    validate(need(length(input$selected_drugs) !=0 , paste0(" Please select < 5 drugs")))  
    
    x <- get_drug_flt_normViab()
    drugViab_dosages <- dcast(x, group+experiment+stage+cellLine+drug+replicate ~ conc, value.var="normViability",
                              fun.aggregate = function(x) mean(x))
    
    m <- drugViab_dosages[,-c(1:6)] * 100  #to convert fraction to percentage viability
    rowAnnotation <- drugViab_dosages[,c('drug'),drop=F]
    #m.scaled <- t(scale(t(m)))
    
    #convert colnames to microMolar
    colnames(m) <- as.numeric(colnames(m))*1e+6
    #cluster_rows = if(drugViability_heatmap_col_cluster == TRUE) TRUE else NA
    aheatmap(m,
             scale='none',
             distfun="euclidean",
             Colv=NA,
             Rowv =  if(input$drugViability_heatmap_col_cluster == TRUE) TRUE else NA,
             annRow = rowAnnotation,
             info=TRUE,
             cexRow=0,
             main = 'Cell Viability v/s Drug Dosage(microMolar)',
             sub = 'color signifies cell viability %'
             )
             
  })
  
  output$drugResponse_plots <- renderPlot({
    
    validate(need(length(input$selected_drugs) != 0, paste0(" Please select drug/s (max upto 4)")))  
    validate(need(length(input$selected_drugs) < 5, paste0(" Please select < 5 drugs")))  
    flt_drug_normViab <- get_drug_flt_normViab()

    doseResp <- ddply(.data=flt_drug_normViab, .variables = c('drug', 'cellLine', 'experiment', 'group'), 
                      .fun = tmp_iterator, .parallel = T)
    
    facet_by <- paste(input$facet_by, collapse = ' + ' )
    facet_by <- formula(paste(facet_by, ' ~ ', input$dose_response_plot_splitBy))
    
    
    color_options = c('drug', 'cellLine')
    color_by <- color_options[!color_options %in% input$dose_response_plot_splitBy]
    
    p <- ggplot(data=doseResp, aes_string(x="fittedX", y="fittedY*100", group=color_by))
    p <- p + geom_line(aes_string(color=color_by)) + facet_grid( facet_by) + theme_bw()
    p <- p + geom_hline(aes(yintercept=50), color='grey50', linetype='dashed')
    p <- p + xlab('micromolar conc') + ylab('cell viability %') 
    p <- p + scale_x_continuous(breaks=seq(from=-10,to=-1,by=1),
                                labels = lapply(seq(from=-10,to=-1,by=1), function(x) (10^x)*(1e+6)  ))
    
    p
  })

})


############################
#Public Data tab code
############################

#   selected_pathways <- reactive({
#     #selected grouped pathways
#     selected_grps <- names(grouped_pathways) %in% input$grouped_pathways
#     pathways <- unlist(grouped_pathways[selected_grps])
#     #selected pathways based on entered genes
#     pathways_based_on_selected_genes <- genes_to_pathways$pathway[genes_to_pathways$genes %in% selected_genes()]
#     #now add any individual pathways
#     pathways <- unique(c(pathways, input$pathways, pathways_based_on_selected_genes))
#     pathways
#   })
#   
#   genes_in_pathways <- reactive({
#     #selected grouped pathways
#     selected_grps <- names(grouped_pathways) %in% input$grouped_pathways
#     pathways <- unlist(grouped_pathways[selected_grps])
#     #now add any individual pathways
#     pathways <- unique(c(pathways, input$pathways))
#     #find the genes involved in all the selected pathways
#     selected_genes_based_on_pathways <- genes_to_pathways[genes_to_pathways$pathway %in% pathways,][['genes']]
#     selected_genes_based_on_pathways
#   })  
#   
#   selected_genes <- reactive({
#     input$searchByGenes_button
#     selected_genes <- isolate(input$custom_gene_list)
#     geneList <- unlist(strsplit(selected_genes,split=c('[\\s+,\\n+\\r+)]'),perl=T))
#     #convert everything to upper case
#     geneList <- toupper(geneList)
#     geneList <- geneList[ !geneList == "" ] #remove the blank entries
#     geneList
#   })
#   
#   selected_samples <- reactive({
#     selected_samples <- NULL
#      if(! is.null(input$selected_geo_studies)){
#        publicData_phenotype_flt <- filter(publicData_phenotype, studyId %in% input$selected_geo_studies)
#        selected_samples <- unique(paste0(publicData_phenotype_flt$studyId, "_", publicData_phenotype_flt$sampleId))
#      }
#    selected_samples
#   })
# 
#   
#   filtered_pathway_enrichment_matrix <- reactive({
#       validate(
#         need(length(selected_pathways()) > 0, "Please select atleast one choice in option 1,2,3")
#       )
#       #filter based on pathways
#       rows_to_keep <- rownames(pathway_enrichment_scores) %in% selected_pathways()
#       filtered_matrix <- pathway_enrichment_scores[rows_to_keep,] 
#       
#       #filter based on selected studies
#       if(! is.null(selected_samples())){
#         rows_to_keep <- colnames(filtered_matrix) %in%  selected_samples()
#         filtered_matrix <- filtered_matrix[,rows_to_keep]
#       }
#       filtered_matrix
#   })
#   
#   
#   filtered_expression_matrix <- reactive({
#     #filter based on genes
#     genes <- as.character(unique(c(selected_genes(), genes_in_pathways())))
#     validate( need(length(genes) > 0, "Please select atleast one choice in option 1,2,3") )
#     rows_to_keep <- rownames(publicData_expMat) %in% genes
#     filtered_expression_matrix <- publicData_expMat[rows_to_keep,]  
#         
#     #filter based on selected studies
#     if(! is.null(selected_samples())){
#       rows_to_keep <- colnames(filtered_expression_matrix) %in%  selected_samples()
#       filtered_expression_matrix <- filtered_expression_matrix[,rows_to_keep]
#     }
# 
#     #filter rows with NA values
#     filtered_expression_matrix <- na.omit(filtered_expression_matrix)
#     
#     #make sure the filtered data only has < 1000 genes
#     validate(  need(nrow(filtered_expression_matrix) < 1000, 
#                paste0("Due to longer compute time Expression clustering will only work upto 1000 
#                        genes \n currently entered ",nrow(filtered_expression_matrix), " genes"))
#                )
#     filtered_expression_matrix  
#   })
#   
#   
#   pubData_filteredAnnotation <- reactive({
#     annotation <- publicData_phenotype[,c('studyId', 'sampleId')]
#     annotation['studyType'] = 'meningioma'
#     annotation$studyType[annotation$studyId %in% get_schwannoma_studies()] = 'schwannoma'
#     
#     #if only one geo study is selected selected 
#     if( length(input$selected_geo_studies) == 1 ){
#       annotation <- filter(publicData_phenotype, studyId %in% input$selected_geo_studies )
#       annotation <- dcast(annotation, sampleId + studyId  ~ variable)
#       cols_to_keep <- c('sampleId', 'studyId', input$pubData_selectedPhenotypes)
#       annotation <- annotation[,cols_to_keep]
#     }  
#     #remove duplicated rows
#     rows_to_keep <- ! duplicated(annotation)
#     annotation <- annotation[rows_to_keep,]
#     row.names(annotation) <- paste0(annotation$studyId, '_', annotation$sampleId)
#     annotation$sampleId <- NULL
#     annotation
#   })
#   
#   get_schwannoma_studies <- reactive({ schwannoma_studies })
#   
#   #heatmap for pubData pathway enrichment 
#   output$pubData_pathway_heatmap <- renderPlot({
#     fontsize_row=10
#     if(nrow(filtered_pathway_enrichment_matrix()) > 100){
#       fontsize_row = 0
#     }
#     scaled_mat = t(scale(t(filtered_pathway_enrichment_matrix())))
#     withProgress(message = "clustering & rendering heatmap, please wait", 
#                  detail = "This may take a few moments...", {
#                             expHeatMap(scaled_mat, pubData_filteredAnnotation(),
#                                        scale=F,
#                                        fontsize_col=0,
#                                        fontsize_row=fontsize_row)
#                 })
#   })
#   
#   #heatmap for EXPRESSION pubData 
#   output$pubData_expression_heatmap <- renderPlot({  
#     fontsize_row=10
#     if(nrow(filtered_expression_matrix()) > 100){
#       fontsize_row = 0
#     }
#     scaled_expmat = t(scale(t(filtered_expression_matrix())))
#      withProgress( message = "clustering & rendering heatmap, please wait", 
#                    detail = "This may take a few moments...", {
#                      expHeatMap(scaled_expmat, 
#                                 pubData_filteredAnnotation(),
#                                 scale = F,
#                                 fontsize_col=0,
#                                 fontsize_row=fontsize_row)               
#                    })
#   })
#  
#   #update the available phenotypes for a single study
#   observe({
#     phenotypes <- filter(publicData_phenotype, studyId %in% input$selected_geo_studies)
#     choices = unique(as.character(phenotypes[['variable']]))
#     updateSelectInput(session = session,
#                       inputId = "pubData_selectedPhenotypes",
#                       choices = choices)
#   })

#######################################################
#END public data tab code
#######################################################



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



#   kinomeData_scatterPlot <- reactive({
#     df <- get_filtered_kinomeData_by_userSelected_thresholds()
#     selected_points <- linked_brush(keys=df$id, "red")
#     
#     vis <- df %>%
#       ggvis(x = ~count, y = ~log2ratio, fill = ~condition, key := ~id, size.hover := 200) %>%
#       layer_points() %>%
#       add_legend(c("fill")) %>%
#       add_tooltip(tooltip_values, "hover") %>%
#       set_options(width=450, height=300) %>%
#       selected_points$input() 
#     
#     return(vis)
#   })
#   kinomeData_scatterPlot %>% bind_shiny("kinomeData_scatterPlot") 




# linked_brush <- function(keys, fill = "red") {
#   stopifnot(is.character(fill), length(fill) == 1)
#   
#   rv <- shiny::reactiveValues(under_brush = character(), keys = character())
#   rv$keys <- isolate(keys)
#   
#   input <- function(vis) {
#     handle_brush(vis, fill = fill, on_move = function(items, ...) {
#       rv$under_brush <- items$key__
#     })
#   }
#   
#   set_keys <- function(keys) { rv$keys <- keys }
#   set_brush <- function(ids) { rv$under_brush <- ids }
#   selected_r <- reactive(rv$keys %in% rv$under_brush)
#   fill_r <- reactive(c("black", fill)[selected_r() + 1])
#   
#   list(
#     input = input,
#     selected = create_broker(selected_r),
#     fill = create_broker(fill_r),
#     set_keys = set_keys,
#     set_brush = set_brush
#   )
# }


#   output$kinome_barPlot <- renderChart({
#     p <- nPlot(ratio ~ Uniprot , data=get_ordered_kinomeData(), group="condition", 
#                type="multiBarChart")
#     p$params$width <- 800
#     p$params$height <- 500
#     p$addParams(dom = 'kinome_barPlot')
#     return(p)
#   })



#get_ordered_kinomeData() %>%  
#   output$down_kinomeBarPlot <- downloadHandler <- (
#     filename = function(){
#       paste('synodos_kinomeBarPlot.png')
#     },
#     content = function(file){
#       png(filename) 
#       
#     }
#     )


#1. user def scatterPlot
#tooltip
# tooltip_values <- function(x){
#   data <- get_filtered_kinomeData_by_userSelected_thresholds()
#   if(is.null(x)) return(NULL)
#   row <- filter(data, id == x$id)
#   selected_cols <- c('Description', 'Gene', 'Uniprot', 'Family',
#                      'ratio', 'variability', 'count', 'uniq_peptides')
#   row <- row[,selected_cols]
#   paste0(names(row), ": ", format(row), collapse = "<br />")
# }




