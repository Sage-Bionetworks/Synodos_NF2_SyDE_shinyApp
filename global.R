#global.R
options(stringsAsFactors = FALSE)
library(synapseClient)
library('rCharts')
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
registerDoMC(4)

synapseLogin()


#source the heatmap code
source_url("https://raw.githubusercontent.com/apratap/apRs/master/expression_heatmap.R")

#loding the global data for publicData 
source("global_publicData.R")

#loading the global kinome screening data
source("global_kinomeScreens.R")

#loading the global DrugScreens Data
source("global_DrugScreens.R")

#source("global_testing.R")



global_cellLines_metadata_link <- 'https://www.synapse.org/#!Synapse:syn2774452/tables/query/%7B%22limit%22:10,%20%22sql%22:%22SELECT%20*%20FROM%20syn2774452%22,%20%22isConsistent%22:true,%20%22offset%22:0%7D'
