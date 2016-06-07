##download the precomputed IC Vals
# drug_ICVals.synId <- 'syn2773794'
# drug_ICVals.file <- synGet(drug_ICVals.synId)
# drug_ICVals <- read.delim(drug_ICVals.file@filePath, check.names=F, sep="\t", header=T)
# drug_ICVals <- filter(drug_ICVals, goodNess_of_fit > .70 & hillSlope < 0 )
# 
# drug_ICVals <- rename(drug_ICVals,sample = cellLine)
# 
# drugScreenData <- drug_ICVals[c("sample","IC50","maxEfficacy","drug")]

# Synodos Data

synodosData <- fread("synodosData.csv")
names(synodosData) <- c("protocolName","drugID","sample","sampleType","nf2Status",
                        "AC50","curveClass","maxResp","logAC50","drug","AUC",
                        "AUCfit","target")

dataCols <- c("sample","AC50","maxResp","curveClass","drugID","drug","target","AUC")
drugScreenData <- as.data.frame(synodosData)[dataCols]

drugScreenData$IC50 <- NA
