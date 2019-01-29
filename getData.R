###########
# Raw Drug Screening Data
############

rawData <- "syn6138251"
rawData <- synGet(rawData)
rawData <- read.delim(rawData$path, check.names=F, sep="\t", header=T)
colnames(rawData)[1] <- "sample"
select_col <- c("sample", "conc", "drug", "replicate", "normViability")
rawData <- rawData[,select_col]

#################
# Summarized Drug Screening Data
#################

summarizedData <- "syn6138237"
summarizedData <- synGet(summarizedData)
summarizedData <- read.delim(summarizedData$path, check.names=F, sep="\t", header=T)
select_col <- c("cellLine","drug","IC50","maxEfficacy","trapezoid")
summarizedData <- summarizedData[,select_col]
colnames(summarizedData) <- c("sample", "drug", "IC50", "maxResp","AUC")

#convert IC50 to uM
summarizedData$IC50 <- summarizedData$IC50*1e+6

#convert maxResp to 1-100 scale
summarizedData$maxResp <- summarizedData$maxResp*100

#add missing columns
summarizedData$AC50 <- NA
summarizedData$target <- NA
summarizedData$curveClass <- NA

#################
# CPM RNASeq Data Matrix (Rank Normalized ExpressionSet)
#################

RNAseq <- "syn10845587"
RNAseq <- synGet(RNAseq)$path
RNAseq <- readRDS(file = RNAseq)

###########
# Pathway List
###########

MSIGDB_syn<-synGet("syn2227979")
load(MSIGDB_syn$path) #available as MSigDB R object

pathways_list <- c(MSigDB$C2.CP.BIOCARTA, MSigDB$C2.CP.KEGG, MSigDB$C2.CP.REACTOME)

#################
# Drug Treated Kinome Ratios Data Matrix
#################

kinometx <- "syn10845736"
kinometx <- readRDS(synGet(kinometx)$path)

#################
# Drug Treated Kinome Ratios Data Table
#################

remove synapse dependency
basekin<-read.table(synGet("syn5840701")$path, sep = "\t", header = TRUE, comment.char = "")


Syn5.Syn1.base <- basekin %>% filter(cellLine=="Syn5", referenceSample=="Syn1") %>% 
  group_by(Gene) %>% 
  dplyr::summarize(mean.log2ratio = mean(log2ratio, na.rm = TRUE), sem = (sd(log2ratio, na.rm = TRUE)/(sqrt(length(log2ratio)))), comp = "Syn5_Syn1")  %>% 
  filter(abs(sem)<abs(mean.log2ratio)) %>% 
  filter(abs(mean.log2ratio)>0.1) 

Syn5.Syn1.base$Gene <- reorder(Syn5.Syn1.base$Gene, Syn5.Syn1.base$mean.log2ratio)

HS01.HS11.base <- basekin %>% filter(cellLine=="HS01", referenceSample=="HS11") %>% 
  group_by(Gene) %>% 
  dplyr::summarize(mean.log2ratio = mean(log2ratio, na.rm = TRUE), sem = (sd(log2ratio, na.rm = TRUE)/(sqrt(length(log2ratio)))), comp = "HS01_HS11") %>% 
  filter(abs(sem)<abs(mean.log2ratio)) %>% 
  filter(abs(mean.log2ratio)>0.1) 

HS01.HS11.base$Gene <- reorder(HS01.HS11.base$Gene, HS01.HS11.base$mean.log2ratio)


