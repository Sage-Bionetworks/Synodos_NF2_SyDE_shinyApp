options(stringsAsFactors = FALSE)

filter_df <- function(df,pattern){
  temp_df <- df[with(df,grepl(pattern,entity.name)),]  
  temp_df['geo_study_id'] = gsub(paste0('_',pattern,'.*'),"",temp_df$entity.name)
  temp_df
}


#get the MsigDB object
cat('Reading the MSIGDB object from synapse...')
MSIGDB_syn<-synGet("syn2227979")
load(MSIGDB_syn@filePath) #available as MSigDB R object
cat('..Done\n\n')
pathways <- c(MSigDB$C2.CP.BIOCARTA, MSigDB$C2.CP.KEGG, MSigDB$C2.CP.REACTOME)


#uniq list of all the genes in the pathways
genes_in_pathways <- unique(as.character(unlist(pathways)))

#data frame to map between genes and pathways
genes_to_pathways <- data.frame('genes'=unlist(pathways) , 'pathway'=names(unlist(pathways)))
#fixing the numbers suffix in the pathways
genes_to_pathways$pathway <- gsub('\\d+','',genes_to_pathways$pathway,perl=T)


############
#download combined pathway enrichment file
############
pathway_file_id <- "syn2580137"
pathway_enrichment_scores <- read.table(synGet(pathway_file_id)@filePath, header=T)
#modify the colnames
temp <- strsplit(colnames(pathway_enrichment_scores),split="_")
new_colnames <- (paste0(gsub('\\..*', '', lapply(temp, "[[", 1)), "_",
                 gsub('\\..*', '', lapply(temp, "[[", 2))
                 ))
colnames(pathway_enrichment_scores) <- new_colnames


##########
#download the combined expression file
##########
exp_file_id <- "syn2626803"
publicData_expMat <- read.table(synGet(exp_file_id)@filePath, header=T)
temp <- strsplit(colnames(publicData_expMat), split="_")
new_colnames <- (paste0(gsub('\\..*', '', lapply(temp, "[[", 1)), "_",
                        gsub('\\..*', '', lapply(temp, "[[", 2))
))
colnames(publicData_expMat) <- new_colnames

#ignore all data from the GSE56157/GSE58037 studies as the values seem to be not properly normalized
#these are the illumina studies
publicData_expMat <- publicData_expMat[,!grepl('GSE56157|GSE58037', colnames(publicData_expMat))]





##########
#download the phenotype file for each study
##########
synQuery_results <- synQuery("select id,name from entity where benefactorId == 'syn2347420'")
publicData_phenotype_files <- filter_df(synQuery_results,'phenotype')

prepare_phenotype_data <- function(df){
  synId <- df[[2]]
  tempDF <- read.table(synGet(synId)@filePath, header=T, sep="\t")
  if( ncol(tempDF) > 1) {
    tempDF <- melt(tempDF, id.vars=c('sampleId'), na.rm=F)
  }
  else{
    tempDF['variable'] = NA
    tempDF['value'] = NA
  }
  tempDF['studyId'] = df[[3]]
  tempDF
}

publicData_phenotype <- apply(publicData_phenotype_files, 1, prepare_phenotype_data)
publicData_phenotype <- Reduce(rbind, publicData_phenotype)


#study types
schwannoma_studies = c('GSE30563','GSE39645')


#################
##Pathways groupings
#################
#pathways with NF2 genes
#pathways_with_NF2_genes <- names(pathways[sapply(pathways, function(x){sum(grepl('NF2',x))}) > 0])
#NF2 related pathways
NF2_related_pathways  <- c( "KEGG_MTOR_SIGNALING_PATHWAY", "KEGG_WNT_SIGNALING_PATHWAY",
                            "REACTOME_SIGNALING_BY_HIPPO" )

#KINASE related pathways
KINASE_related_pathways <- names(pathways)[grep('KINASE',names(pathways))]


grouped_pathways <- list('NF2_related_pathways' = NF2_related_pathways,
                         'KINASE_related_pathways' = KINASE_related_pathways)



