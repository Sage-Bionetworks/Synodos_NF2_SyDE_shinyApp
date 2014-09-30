#global file for Kinome screening data 

test_MIB_MS_datasheet <- "syn2506526"
test_MIB_MS_datasheet <- synGet(test_MIB_MS_datasheet)@filePath
df <- read.xls(test_MIB_MS_datasheet, sheet=2,header=T,stringsAsFactors=FALSE)
df <- df[,c(1:20)]


#get the colnames which are other than the first 5 cols and not containing the variability and count
ratio_cols_match_pattern = paste0("Variability|Count|", paste(colnames(df)[1:5], collapse="|"))
ratio_cols <- colnames(df)[!grepl(ratio_cols_match_pattern,colnames(df),perl=T)]
variability_cols <- colnames(df)[grepl("Variability",colnames(df),perl=T)]
count_cols <- colnames(df)[grepl("Count",colnames(df),perl=T)]

ratios <- melt(df, id.vars=colnames(df)[1:5], measure.vars = ratio_cols, variable.name = 'sample' , value.name = 'ratio')
ratio_variability <- melt(df,id.vars=colnames(df)[1:5], measure.vars = variability_cols, variable.name = 'sample' , value.name = 'variability')
ratio_variability$sample <- gsub('\\.Var.*','', ratio_variability$sample,perl=T)
psm_counts <-   melt(df,id.vars=colnames(df)[1:5], measure.vars = count_cols, variable.name = 'sample' , value.name = 'count')
psm_counts$sample <- gsub('\\.Count.*','', psm_counts$sample,perl=T)

#create the final data frame
merge.all <- function(x,y){ merge(x,y) }
kinomeData <-Reduce(merge.all, list(ratios, ratio_variability, psm_counts))


filtered_data <- kinomeData[!apply(kinomeData, 1, function(x) {sum(is.na(x)) > 0}),]
filtered_data['log2_ratio'] = log2(filtered_data$ratio)

head(filtered_data)
