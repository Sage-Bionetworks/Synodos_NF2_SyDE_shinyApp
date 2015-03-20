#global file for Kinome screening data 


temp_kinomeData_dataProcessing <- function(df, iTRAQ_to_cellLine,runCondition ){
  
  #get the colnames which are other than the first 5 cols and not containing the variability and count
  ratio_cols_match_pattern = paste0("Variability|Count|", paste(colnames(df)[1:5], collapse="|"))
  ratio_cols <- colnames(df)[!grepl(ratio_cols_match_pattern,colnames(df),perl=T)]
  variability_cols <- colnames(df)[grepl("Variability",colnames(df),perl=T)]
  count_cols <- colnames(df)[grepl("Count",colnames(df),perl=T)]
  
  ratios <- melt(df, id.vars=colnames(df)[1:5], measure.vars = ratio_cols, variable.name = 'sample' , value.name = 'ratio')
  ratio_variability <- melt(df,id.vars=colnames(df)[1:5], measure.vars = variability_cols, variable.name = 'sample' , value.name = 'variability')
  ratio_variability$sample <- gsub(' Var.*','', ratio_variability$sample,perl=T)
  psm_counts <-   melt(df,id.vars=colnames(df)[1:5], measure.vars = count_cols, variable.name = 'sample' , value.name = 'count')
  psm_counts$sample <- gsub(' Count.*','', psm_counts$sample,perl=T)
  
  #create the final data frame
  merge.all <- function(x,y){ merge(x,y) }
  kinomeData <-Reduce(merge.all, list(ratios, ratio_variability, psm_counts))
  condition <- lapply(strsplit(as.character(kinomeData$sample), split='/'),
                      function(x) paste0(iTRAQ_to_cellLine[[x[[1]]]], '/', iTRAQ_to_cellLine[[x[[2]]]]) )
  condition <- paste0(condition, ' ', runCondition)
  kinomeData['condition'] <- condition
  kinomeData
}



#process KINOME RUN 1 : FullSerum
kinomeRun1 <- synGet('syn2679231')
kinomeRun1 <- read.xls(kinomeRun1@filePath, sheet=2, header=T, check.names=F)
#keep only selected columns (75% co-isolation)
kinomeRun1 <- kinomeRun1[,c(1:14)]
kinomeRun1_iTRAQ_to_cellLine <- list('116'='Syn2','117'='Syn3','114'='Syn5_1','115'='Syn5_2')                                                    
kinomeRun1_condition <- 'FullSerum'
kinomeRun1 <- temp_kinomeData_dataProcessing(kinomeRun1, kinomeRun1_iTRAQ_to_cellLine, kinomeRun1_condition)


#process KINOME RUN 2 : Serum Free
kinomeRun2 <- synGet('syn2679230')
kinomeRun2 <- read.xls(kinomeRun2@filePath, sheet=2, header=T, check.names=F)
#keep only selected columns (75% co-isolation)
kinomeRun2 <- kinomeRun2[,c(1:14)]
kinomeRun2_iTRAQ_to_cellLine <- list('114'='Syn2','115'='Syn3','116'='Syn5_1','117'='Syn5_2')                                                    
kinomeRun2_condition <- 'SerumFree'
kinomeRun2 <- temp_kinomeData_dataProcessing(kinomeRun2, kinomeRun2_iTRAQ_to_cellLine, kinomeRun2_condition)

#combine all the kinome runs
kinomeData <- rbind(kinomeRun1, kinomeRun2)

kinomeData$variability <- as.numeric(kinomeData$variability)
kinomeData['log2ratio'] = log2(kinomeData$ratio)
kinomeData['log2ratio_max'] =  kinomeData$log2ratio    + kinomeData$log2ratio*(kinomeData$variability/100)
kinomeData['log2ratio_min'] =  kinomeData$log2ratio    - kinomeData$log2ratio*(kinomeData$variability/100)
kinomeData['uniq_peptides'] = kinomeData['# Unique Peptides']

