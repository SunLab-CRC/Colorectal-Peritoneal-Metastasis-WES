<PROJECT_ROOT>
<PROJECT_ROOT>
setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
rm(list=ls())
amp <- read.csv('TCGA_2012Nature_SCNA_Amp.csv')
sum(amp$GenesinPeak)
scna_amp <- c()
for(i in 1:nrow(amp)){
	sgen <- strsplit(amp$Candidates[i],', ')[[1]]
	scna_amp <- c(scna_amp,sgen)
	rm(sgen, i)
}
scna_amp <- gsub('\\[','',scna_amp)
scna_amp <- gsub('\\]','',scna_amp)
scna_amp <- scna_amp[!duplicated(scna_amp)]
write.table(scna_amp, file='TCGA_2012Nature_SCNA_Amp_candidate_Gene.txt',row.names=F,quote=F,col.names=F,sep='\t')
rm(list=ls())
amp <- read.csv('TCGA_2012Nature_SCNA_Del.csv')
amp <- subset(amp, !is.na(genesinpeak))
sum(amp$genesinpeak)
scna_amp <- c()
for(i in 1:nrow(amp)){
	sgen <- strsplit(amp$Candidates[i],', ')[[1]]
	scna_amp <- c(scna_amp,sgen)
	rm(sgen, i)
}
scna_amp <- gsub('\\[','',scna_amp)
scna_amp <- gsub('\\]','',scna_amp)
scna_amp <- scna_amp[!duplicated(scna_amp)]
write.table(scna_amp, file='TCGA_2012Nature_SCNA_Del_candidate_Gene.txt',row.names=F,quote=F,col.names=F,sep='\t')
setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
rm(list=ls())
rm(list=ls())
cgc <- read.csv('COSMIC_CancerGeneCensus_v97/COSMIC_SCNA_driver_Amp.csv')
cgc <- cgc[,c('Gene.Symbol','Synonyms.1')]
nat <- read.delim('2012_Nature_TCGA_CRC/TCGA_2012Nature_SCNA_Amp_candidate_Gene.txt',h=F)
colnames(nat) <- 'Symbol'
nat$ID <- NA
cgc$database <- 'COSMIC'
nat$database <- 'TCGA_Nature2012'
amp <- rbind(cgc,nat)
write.csv(amp, file = "SCNA_driver_CRC_Amp_candidate_Gene.csv", row.names=F,quote=F)
rm(list=ls())
cgc <- read.csv('COSMIC_CancerGeneCensus_v97/COSMIC_SCNA_driver_Del.csv')
cgc <- cgc[,c('Gene.Symbol','Synonyms.1')]
nat <- read.delim('2012_Nature_TCGA_CRC/TCGA_2012Nature_SCNA_Del_candidate_Gene.txt',h=F)
colnames(nat) <- 'Symbol'
nat$ID <- NA
cgc$database <- 'COSMIC'
nat$database <- 'TCGA_Nature2012'
amp <- rbind(cgc,nat)
write.csv(amp, file = "SCNA_driver_CRC_Del_candidate_Gene.csv", row.names=F,quote=F)
