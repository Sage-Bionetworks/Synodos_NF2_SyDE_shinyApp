#global file for Kinome screening data 


#get the global kinome data
kinomeData_synid <- 'syn4259365'
kinomeData <- synGet(kinomeData_synid)
kinomeData <- read.table(kinomeData@filePath, sep="\t", header=T)  


