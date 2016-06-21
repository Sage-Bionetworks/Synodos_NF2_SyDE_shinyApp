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
library("nplr")
#require(memoise)
library("devtools")
library("ggplot2")
library("data.table")
library("doMC")
library("NMF")
library("gridExtra")
library("futile.logger")
registerDoMC(4)

flog.threshold(DEBUG, name='server')
flog.threshold(DEBUG, name='ui')
flog.threshold(DEBUG, name='global')
flog.threshold(INFO, name='synapse')

synapseLogin()

flog.debug("Starting App...", name="server")

flog.debug("Loading module...", name="server")
#source("drugScreenModule.R")
source_https <- function(url, ...) {
  # load package
  require(RCurl)
 
  # parse and evaluate each .R script
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}

source_https("https://raw.githubusercontent.com/Sage-Bionetworks/shinyModules/master/drugScreen/drugScreenModule.R?token=APzGNGtNcc2N3bFPsW74vTsq-6DFrEilks5XaufuwA%3D%3D")

flog.debug("Loading data...", name="server")
source("getData.R")
