############
# Raw Data
############
#download the DMSO norm data for MGH and UCF Drug Screens
UCF_normViab <- 'syn2773870'
UCF_normViab <- synGet(UCF_normViab)
UCF_normViab <- read.delim(UCF_normViab@filePath, check.names=F, sep="\t", header=T)
UCF_normViab$cellLine <- gsub("^ ", "", UCF_normViab$cellLine)
UCF_normViab$cellLine <- gsub("Nf2 --", "Nf2--", UCF_normViab$cellLine)

MGH_normViab <- 'syn2773792'
MGH_normViab <- synGet(MGH_normViab)
MGH_normViab <- read.delim(MGH_normViab@filePath, check.names=F, sep="\t", header=T)
MGH_normViab$cellLine <- gsub("^ ", "", MGH_normViab$cellLine)

#drop unnecassary cols
drop_cols <- c('plate', 'medianDMSO', 'viability', 'experiment','stage')
MGH_normViab <- MGH_normViab[, !colnames(MGH_normViab) %in% drop_cols]
UCF_normViab <- UCF_normViab[, !colnames(UCF_normViab) %in% drop_cols]

#align the columns of the two dataframe before combining them
UCF_normViab <- UCF_normViab[,colnames(MGH_normViab)]

#combined drug normViab
#rbindlist is fastest for concatenating data.frames/ data.tables by row
rawData <- as.data.frame(rbindlist(list(UCF_normViab, MGH_normViab))) 
colnames(rawData)[1] <- "sample"


#################
# Summarized Data
#################

#download the precomputed IC Vals for MGH and UCF Data
UCF_ICvals <- 'syn2773891'
UCF_ICvals <- synGet(UCF_ICvals)
UCF_ICvals <- read.delim(UCF_ICvals@filePath, check.names=F, sep="\t", header=T)

MGH_ICvals <- 'syn2773794'
MGH_ICvals <- synGet(MGH_ICvals)
MGH_ICvals <- read.delim(MGH_ICvals@filePath, check.names=F, sep="\t", header=T)

#align the columns of the two dataframe before combining them
MGH_ICvals <- UCF_ICvals[,colnames(UCF_ICvals)]

#combined 
summarizedData <- as.data.frame(rbindlist(list(MGH_ICvals, UCF_ICvals)))
select_col <- c("cellLine","drug","IC50","maxEfficacy")
summarizedData <- summarizedData[,select_col]
colnames(summarizedData) <- c("sample", "drug", "IC50", "maxResp")
summarizedData$maxResp <- summarizedData$maxResp*100
summarizedData$AC50 <- NA
summarizedData$target <- NA
summarizedData$AUC <- NA
summarizedData$curveClass <- NA

