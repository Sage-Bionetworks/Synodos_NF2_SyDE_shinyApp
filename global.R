#global.R
options(stringsAsFactors = FALSE)
library(stringr)
library(synapseClient)
library("shiny")
library("rbokeh")
library("htmlwidgets")
library("RCurl")
library("reshape2")
library("scales")
library("gdata")
library("plyr")
library("dplyr")
library("devtools")
library("ggplot2")
library("data.table")
library("doMC")
library("NMF")
library("shinyBS")
library("ggvis")
registerDoMC(4)
synapseLogin()

#source the heatmap code
source_url("https://raw.githubusercontent.com/apratap/apRs/master/expression_heatmap.R")

#loding the global data for publicData 
#cat('loading global public data ...')
#source("global_publicData.R")
#cat('Done \n\n')

cat('loading the MSigDB pathway <-> genes map....')
MSIGDB_syn<-synGet("syn2227979")
load(MSIGDB_syn@filePath) #available as MSigDB R object
#pathways_list <- c(MSigDB$C2.CP.BIOCARTA, MSigDB$C2.CP.KEGG, MSigDB$C2.CP.REACTOME)
pathways_list <- c(MSigDB$C2.CP.KEGG)
global_pathway_gene_map <- ldply(pathways_list, function(x) {  
                                            data.frame(gene=as.character(x))})
colnames(global_pathway_gene_map) <- c('pathway', 'gene')
cat('Done \n\n')

#loading the global kinome screening data
cat('loading kinome screening data ...')
source("global_kinomeScreens.R")
cat("Done \n\n")

#loading the global DrugScreens Data
cat('loading global drug Screen data ...')
source("global_DrugScreens.R")
cat("Done \n\n")


#source("global_testing.R")
global_cellLines_metadata_link <- 'https://www.synapse.org/#!Synapse:syn3275123/tables/query/eyJsaW1pdCI6MjUsICJzcWwiOiJTRUxFQ1QgKiBGUk9NIHN5bjMyNzUxMjMiLCAiaXNDb25zaXN0ZW50Ijp0cnVlLCAib2Zmc2V0IjowfQ=='
global_drug_metadata_link <- 'https://www.synapse.org/#!Synapse:syn3105963/tables/query/eyJsaW1pdCI6MjUsICJzcWwiOiJTRUxFQ1QgKiBGUk9NIHN5bjMxMDU5NjMiLCAiaXNDb25zaXN0ZW50Ijp0cnVlLCAib2Zmc2V0IjowfQ=='