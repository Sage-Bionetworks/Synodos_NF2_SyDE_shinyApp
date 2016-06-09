drugScreenModule <- function(input,output,session,summarizedData = NULL, rawData = NULL){
  plot_names <- QC_plot_name(summarizedData)
  ns <- NS(id)
  
  output$target_class <- renderUI({
    if(!all(is.na(summarizedData$target))){
      tagList(
        h5("Select by target class"),
        selectInput(ns('selected_class'),NULL, choices = c('ALL', unique(summarizedData$target)),
                    selectize=T, multiple=T)
      )
    }
  })
  
  output$filter1 <- renderUI({
    if(!all(is.na(summarizedData$AC50))){
      tagList(
        sliderInput(ns('ac50_filter'), 'AC50', 
                    min=10, max=100, value=c(1, 1000), 
                    step=5, round=0)
      )
    }
  })
  
  output$filter2 <- renderUI({
    if(!all(is.na(summarizedData$curveClass))){
      tagList(
        sliderInput(ns('curveClass'), 'Curve Class', 
                    min=-10, max=10, value=sort(unique(drugScreenData$curveClass)), 
                    step=1, round=0)
      )
    }
  })
  
  get_selected_samples <- reactive({
    samples <- input$samples
    samples
  })
  
  get_selected_drugs <- reactive({
    drugNames <- input$selected_drugs
    targetClass <- if(!all(is.na(summarizedData$target))) unique(summarizedData[summarizedData$target ==input$selected_class,]$drug) else NA
    drugs <- union(drugNames,targetClass)
    drugs <- drugs[!is.na(drugs)]
    validate(need(length(drugs) != 0, "At least one drug or target class needs to be selected"))
    drugs
  })
  
  get_drug_data <- eventReactive(input$updateButton,{
    validate(need(!is.null(input$samples), "At least one sample needs to be selected" ))
    validate(need(length(input$samples) <= 5, "You can select up to 5 samples" ))
    validate(need(!is.null(input$selected_drugs), "At least one sample needs to be selected" ))
    flt_drug_data <- filter(summarizedData, drug %in% get_selected_drugs())  
    flt_drug_data <- filter(flt_drug_data, sample %in% get_selected_samples())  
    return(flt_drug_data)
  })
  
  # update filter values
  observe({
    drug_data <- get_drug_data()
    
    # max response
    # if(!all(is.na(summarizedData$maxResp))){
    mR_min <- floor(min(drug_data$maxResp,na.rm = T))
    mR_max <- ceiling(max(drug_data$maxResp,na.rm = T))
    updateSliderInput(session, "maxR_filter", min = mR_min, max = mR_max, value = c(mR_min,mR_max))
    #  }
    # AC50 slider
    #     if(!all(is.na(summarizedData$AC50))){
    #       ac50_min <- floor(min(drug_data$AC50,na.rm = T))
    #       ac50_max <- ceiling(max(drug_data$AC50,na.rm = T))
    #       updateSliderInput(session, "ac50_filter", max = ac50_max, value = c(ac50_min,ac50_max))
    #     }
    # IC50 slider
    #if(!all(is.na(summarizedData$IC50))){
    ic50 <- drug_data$IC50
    ic50 <- ic50[!is.na(ic50)]
    ic50 <- ic50[!is.infinite(ic50)]
    ic50_min <- floor(min(ic50))
    ic50_max <- ceiling(max(ic50))
    updateSliderInput(session, "ic50_filter", min = ic50_min, max = ic50_max, 
                      value = c(ic50_min, ic50_max), step = floor((ic50_max - ic50_min)/5))
    #}
    #     # curve class
    #      if(!all(is.na(summarizedData$curveClass))){
    #       cc_min <- floor(min(drug_data$curveClass,na.rm = T))
    #       cc_max <- ceiling(max(drug_data$curveClass,na.rm = T))
    #       updateSliderInput(session, "curveClass", min = cc_min, max = cc_max, value = c(cc_min,cc_max))
    #      }
  })
  
  get_filtered_drug_data <- reactive({
    drug_data <- get_drug_data()
    filtered_data <- drug_data[drug_data$maxResp >= input$maxR_filter[1] & drug_data$maxResp <= input$maxR_filter[2],]
    #if(!all(is.na(summarizedData$curveClass))){
    #  filtered_data <- drug_data[drug_data$curveClass >= input$curveClass[1] & drug_data$curveClass <= input$curveClass[2],]
    #}
    #if(!all(is.na(summarizedData$AC50))){
    #  filtered_data <- drug_data[drug_data$AC50 >= input$ac50_filter[1] & drug_data$AC50 <= input$ac50_filter[2],]
    # }
    # if(!all(is.na(summarizedData$IC50))){
    filtered_data <- filtered_data[filtered_data$IC50 >= input$ic50_filter[1] & filtered_data$IC50 <= input$ic50_filter[2],]
    #  }
    # filtered_data <- filtered_data[!is.na(filtered_data$sample),]
    filtered_data
  })
  
  x_angle <- reactive({
    flt_drug_data <- get_filtered_drug_data()
    num_drugs <- length(unique(flt_drug_data$drug))
    if(num_drugs < 6){
      return(c(0,0.5))
    }else if(num_drugs < 31){
      return(c(30+num_drugs,1))
    }
    return(c(90,1))
  })
  
  # Max Response Plot
  output$drug_max_resp <- renderPlot({
    flt_drug_data <- get_filtered_drug_data()
    drug_levels <- flt_drug_data %>%
      group_by(drug) %>%
      summarise(med=median(maxResp, na.rm=T)) %>%
      arrange(desc(med)) %>% select(drug)
    drug_levels <- drug_levels$drug
    flt_drug_data$drug <- factor(flt_drug_data$drug,levels=drug_levels)
    p <- ggplot(data=flt_drug_data, aes(x=drug, y=maxResp, group=sample)) 
    p <- p + geom_point(aes(color=sample), size=3) + theme_bw(base_size = 15)
    p <- p + theme(text = element_text(size=20), axis.text.x=element_text(angle=x_angle()[1], hjust=x_angle()[2])) + xlab('Drug') + ylab('Response')
    p
  })
  
  # IC50 plot
  output$drugScreen_IC50_plot <- renderPlot({
    flt_drug_data <- get_filtered_drug_data()
    #remove NA and Inf
    flt_drug_data <- flt_drug_data[! is.na(flt_drug_data$IC50), ]
    flt_drug_data <- flt_drug_data[! is.infinite(flt_drug_data$IC50), ]
    #convert to log10
    flt_drug_data$IC50 <- log10(as.numeric(flt_drug_data$IC50)) 
    drug_levels <- flt_drug_data %>%
      group_by(drug) %>%
      summarise(med=median(IC50, na.rm=T)) %>%
      arrange(desc(med)) %>% select(drug)
    drug_levels <- drug_levels$drug
    flt_drug_data$drug <- factor(flt_drug_data$drug,levels=drug_levels)
    p <- ggplot(data=flt_drug_data, aes(x=drug, y=IC50, group=sample)) 
    p <- p + geom_point(aes(color=sample), size=3) + theme_bw(base_size = 15)
    p + theme(text = element_text(size=20), axis.text.x=element_text(angle=x_angle()[1], hjust=x_angle()[2])) + xlab('Drug') + ylab('log 10 IC50')
  })
  
  # AC50 plot
  output$drugScreen_AC50_plot <- renderPlot({
    validate(need("AC50" %in% plot_names, "AC50 data does not exist." ))
    flt_drug_data <- get_filtered_drug_data()
    #remove NA
    flt_drug_data <- flt_drug_data[! is.na(flt_drug_data["AC50"]), ]
    #convert to log10
    flt_drug_data["AC50"] <- log10(as.numeric(flt_drug_data[,"AC50"])) 
    
    drug_levels <- flt_drug_data %>%
      group_by(drug) %>%
      summarise(med=median("AC50", na.rm=T)) %>%
      arrange(desc(med)) %>% select(drug)
    drug_levels <- drug_levels$drug
    flt_drug_data$drug <- factor(flt_drug_data$drug,levels=drug_levels)
    p <- ggplot(data=flt_drug_data, aes_string(x="drug", y="AC50", group="sample")) 
    p <- p + geom_point(aes(color=sample), size=3) + theme_bw(base_size = 15)
    p + theme(text = element_text(size=20), axis.text.x=element_text(angle=x_angle()[1], hjust=x_angle()[2])) + xlab('Drug') + ylab('AC50')
  })
  
  # QC plots
  QC_plot_list <- reactive({
    plotlist <- list()
    flt_drug_data <- get_filtered_drug_data()
    plotlist <- lapply(plot_names, function(x){
      p <- ggplot(data=flt_drug_data, aes_string(x=x, fill="sample")) + geom_density(alpha=.7)
      p + theme_bw(base_size = 15) 
    })
    
    do.call(grid.arrange, c(plotlist, list(ncol = 2)))
  })
  
  output$QC_plots <- renderPlot({
    QC_plot_list()
  })
  
  
  # Dose Response plot    
  output$doseResp_plot <- renderPlot({
    drug_data <- get_filtered_drug_data()
    drugs <- unique(drug_data$drug)
    if(length(drugs) > 8){
      drugs <- drugs[1:8]
    }
    samples <- unique(drug_data$sample)
    normData <- rawData[rawData$sample %in% samples & rawData$drug %in% drugs,]
    var <- c('drug', 'sample')
    #     if(input$replicate){
    #       var <- c(var,"replicate")
    #     }
    doseRespData <- ddply(.data=normData, .variables = var,.fun = tmp_iterator, .parallel = T)
    doseRespData$grp <- doseRespData$sample
    #     if(input$replicate){
    #       doseRespData$grp <- paste(doseRespData$sample,doseRespData$replicate)
    #     }
    p <- ggplot(normData, aes(x = log10(conc*(1e+6)), y = normViability*100)) 
    p <- p + geom_point(aes_string(color="sample")) 
    p <- p + scale_color_brewer(type = "qual", palette = 2, direction = 1)
    p <- p + geom_line(data = doseRespData, aes(x = fittedX+6, y = fittedY*100, colour = sample, group = grp))
    p <- p + facet_wrap(~ drug, ncol = 4) + theme_bw(base_size = 15)
    p <- p + geom_hline(aes(yintercept=50), color='grey50', linetype='dashed')
    p <- p + xlab('log10 micromolar conc') + ylab('cell viability %') 
    p
  })
  
  # Data table
  output$drugScreen_dataTable <- renderDataTable(
    get_filtered_drug_data()
  )
  
  output$downloadData <- downloadHandler(
    filename = function() { 'summarizedData.csv' },
    content = function(file) {
      write.csv(get_filtered_drug_data(), file)
    }
  )
}

shinyServer(function(input, output, session) {
  callModule(drugScreenModule,id = "demo",session = session, summarizedData = summarizedData, rawData = rawData)
})
