############
# Raw Data
############

rawData <- "syn6138251"
rawData <- synGet(rawData)
rawData <- read.delim(rawData@filePath, check.names=F, sep="\t", header=T)
colnames(rawData)[1] <- "sample"
select_col <- c("sample", "conc", "drug", "replicate", "normViability")
rawData <- rawData[,select_col]

#################
# Summarized Data
#################

summarizedData <- "syn6138237"
summarizedData <- synGet(summarizedData)
summarizedData <- read.delim(summarizedData@filePath, check.names=F, sep="\t", header=T)

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

