#ShinyApp : Synodos Data Explorer
#Tab Panel 2 : exploring kinome screening data

library(synapseClient)
library(dplyr)
synapseLogin()
library(reshape2)
library("gdata")


test_MIB_MS_datasheet <- "syn2506526"
test_MIB_MS_datasheet <- synGet(test_MIB_MS_datasheet)@filePath
df <- read.xls(test_MIB_MS_datasheet, sheet=2,header=T)
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
finalData <-Reduce(merge.all, list(ratios, ratio_variability, psm_counts))


group_by(finalData, "Gene") %>% summarise(max(ratio))

  finalData %.% 
  group_by("Gene") %.%
  summarise(count=n())

sessionInfo()
  
  summarise(ratio_sum = sum(log2(ratio)))



                        arrange(desc(ratio_sum))



matches <- match(finalData$Gene, ratioSum_perGene$Gene )
new_order <- c()
for(i in 1:length(matches)){
  new_order <- c(new_order,which(matches == i))
}

finalData <- finalData[new_order,]

missing_data <- finalData[apply(finalData, 1, function(x) {sum(is.na(x)) > 0}),]
filtered_data <- finalData[!apply(finalData, 1, function(x) {sum(is.na(x)) > 0}),]
filtered_data['log2_ratio'] = log2(filtered_data$ratio)




library(ggplot2)
ggplot(filtered_data, aes(x=log(count), fill=sample)) + geom_histogram(position="identity")





library(ggvis)
colnames(finalData)

finalData %>%
  ggvis(x=~count, y=~variability) %>%
  layer_points(fill = ~factor(sample))


finalData %>%
  ggvis(~log2(ratio[!is.na(ratio)])) %>%
  layer_histograms(binwidth = input_slider(1,100,step=1))





filtered_data %>%
  ggvis(y=~log2(ratio), x=~Uniprot) %>%
  layer_bars()


library('rCharts')
dim(filtered_data)
colnames(filtered_data)
nPlot(log2_ratio ~ Uniprot , data=filtered_data, group="sample", type="multiBarChart")
p$params$facet = "family"
p$params$setLib(lib="libraries")
rPlot(log2_ratio ~ Uniprot, data=filtered_data[1:400,], type="bar", color="sample")


colnames(filtered_data)


nPlot(count, data=filtered_data, group="sample", type="stackedAreaChart")








rPlot
filtered_data$Family <- as.character(filtered_data$Family)
table(filtered_data$Family)
