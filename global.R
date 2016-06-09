#global.R
options(stringsAsFactors = FALSE)
library(shiny)
library(shinyIncubator)
library(shinydashboard)
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
library("NMF")
library("gridExtra")
registerDoMC(4)

synapseLogin()

source("memoised_drugScreen.R")
#source("testModule.R")
#source("drugScreenModule.R")

source("getData.R")

QC_plot_name <- function(data){
  result <- "maxResp"
  if(!all(is.na(data$AUC))){
    result <- c(result, "AUC")
  }
  if(!all(is.na(data$AC50))){
    result <- c(result, "AC50")
  }
  return(result)
}