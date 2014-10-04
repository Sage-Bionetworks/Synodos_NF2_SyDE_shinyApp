Drug_normViab <- synGet("syn2753291")
Drug_normViab <- read.delim(Drug_normViab@filePath, check.names=F, sep="\t", header=T)


Drug_ICVals <- synGet("syn2753295")
Drug_ICVals <- read.delim(Drug_ICVals@filePath, check.names=F, sep="\t", header=T)
