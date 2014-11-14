library("nplr")

get_drugResponse_stats <- function(conc,viability,...){
  res <- nplr(conc, viability,...)
  results <- getAUC(res)
  results['goodNess_of_fit'] <- getGoodness(res)
  results['stdErr'] <- getStdErr(res)
  ICx_est = getEstimates(res, targets= c(.10,.20,.30,.40,.50,.60,.70,.80,.90))
  results['IC10'] = ICx_est[1,'x']
  results['IC20'] = ICx_est[2,'x']
  results['IC30'] = ICx_est[3,'x']
  results['IC40'] = ICx_est[4,'x']
  results['IC50'] = ICx_est[5,'x']
  results['IC60'] = ICx_est[6,'x']
  results['IC70'] = ICx_est[7,'x']
  results['IC80'] = ICx_est[8,'x']
  results['IC90'] = ICx_est[9,'x']
  results['bottom_asymptote'] = res@pars['bottom']
  results['top_asymptote'] = res@pars['top']
  results['hillSlope'] =  res@pars['scal']
  
  fittedVals <- data.frame(fittedX = getXcurve(res),
                           fittedY = getYcurve(res))
  results <- cbind(results,fittedVals)
  results
}

tmp_iterator <- function(df){
  tryCatch({
    stats <- get_drugResponse_stats(df$conc, df$normViability, useLog=T)  
  },error=function(e){
    print(dim(df))
    print(df$conc)
    print(df$normViability)
    print(unique(df$cellLine))
    print(unique(df$drug))
    print(unique(df$experiment))
    stop('stopped')
  })
}


#download the DMSO norm data for MGH and UCF Drug Screens
UCF_normViab <- 'syn2773870'
UCF_normViab <- synGet(UCF_normViab)
UCF_normViab <- read.delim(UCF_normViab@filePath, check.names=F, sep="\t", header=T)
UCF_normViab['group'] <- 'UCF'

MGH_normViab <- 'syn2773792'
MGH_normViab <- synGet(MGH_normViab)
MGH_normViab <- read.delim(MGH_normViab@filePath, check.names=F, sep="\t", header=T)
MGH_normViab['group'] <- 'MGH'



#drop unnecassary cols
drop_cols <- c('plate', 'meanDMSO', 'viability')
MGH_normViab <- MGH_normViab[, !colnames(MGH_normViab) %in% drop_cols]
UCF_normViab <- UCF_normViab[, !colnames(UCF_normViab) %in% drop_cols]

#align the columns of the two dataframe before combining them
UCF_normViab <- UCF_normViab[,colnames(MGH_normViab)]

#combined drug normViab
drug_normViab <- as.data.frame(rbindlist(list(UCF_normViab, MGH_normViab))) #rbindlist is fastest for concatenating data.frames/ data.tables by row


#download the precomputed IC Vals for MGH and UCF Data
UCF_ICvals <- 'syn2773891'
UCF_ICvals <- synGet(UCF_ICvals)
UCF_ICvals <- read.delim(UCF_ICvals@filePath, check.names=F, sep="\t", header=T)
UCF_ICvals['group'] = 'UCF'

MGH_ICvals <- 'syn2773794'
MGH_ICvals <- synGet(MGH_ICvals)
MGH_ICvals <- read.delim(MGH_ICvals@filePath, check.names=F, sep="\t", header=T)
MGH_ICvals['group'] = 'MGH'

#align the columns of the two dataframe before combining them
UCF_ICvals <- UCF_ICvals[,colnames(MGH_ICvals)]

#combined 
drug_ICVals <- as.data.frame(rbindlist(list(MGH_ICvals, UCF_ICvals))) #rbindlist is fastest for concatenating data.frames/ data.tables by row
drug_ICVals <- filter(drug_ICVals, goodNess_of_fit > .70 & hillSlope < 0 )



# selected_IC_value <- 50
# ICx <- eval(paste0('IC', selected_IC_value))
# 
# drug_ICVals[ICx] <- log10(as.numeric(drug_ICVals[,ICx]))
# drug_levels <- (drug_ICVals %>%
#                               group_by(drug) %>%
#                               summarise(mean=mean(IC50, na.rm=T)) %>%
#                               arrange(desc(mean)) %>% select(drug))
# drug_levels <- drug_levels$drug
# drug_ICVals$drug <- factor(drug_ICVals$drug,levels=drug_levels)
# p <- ggplot(data=drug_ICVals, aes_string(x="drug", y="IC50", group="cellLine")) + geom_line(aes(color=cellLine)) 
# p <- p + facet_grid(  ~ ) + geom_point(color='grey50') + theme_bw()
# p + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab('Drug') + ylab(paste0(ICx, ' (log10 molar conc)'))
# 

