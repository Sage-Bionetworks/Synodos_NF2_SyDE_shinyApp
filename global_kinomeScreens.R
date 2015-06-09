#global file for Kinome screening data 


#get the global kinome data
kinomeData_synid <- 'syn4259365'
kinomeData <- synGet(kinomeData_synid)
kinomeData <- read.table(kinomeData@filePath, sep="\t", header=T)  


#get the proteins pval
kinome_proteins_pvals_synid <- 'syn4328906'
kinome_proteins_pvals <- synGet(kinome_proteins_pvals_synid)
kinome_proteins_pvals <- read.table(kinome_proteins_pvals@filePath, sep="\t", header=T)  

