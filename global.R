#global.R
options(stringsAsFactors = FALSE)
library(synapseClient)
library('rCharts')
library("RCurl")
library("reshape2")
library("scales")
library("gdata")
library("dplyr")
synapseLogin()



#load the heatmap code
source("heatmap_code.R")

#loding the global data for publicData 
source("global_publicData.R")

#loading the global kinome screening data
source("global_kinomeScreens.R")

#source("global_testing.R")


get_synodos_banner <- function(...){
  img(src="synodos-banner.jpg",...)
}

 
# library("org.Hs.eg.db")
# columns(org.Hs.eg.db)
# k <- keys(org.Hs.eg.db,keytype="SYMBOL")
# 
# 
# hg19_gene_annot <- select(org.Hs.eg.db, keys=k, columns=c("GENENAME","ALIAS", "CHR", "CHRLOC"), keytype="SYMBOL")
# head(hg19_gene_annot)
# x <- hg19_gene_annot %>%
#   group_by(SYMBOL, GENENAME, CHR, CHRLOC) %>%
#   summarise(paste(ALIAS, collapse=', ')) %>%
#   head(5)
