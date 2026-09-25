options(stringsAsFactors=F)
library(data.table)
setwd('<PROJECT_ROOT>')
rm(list=ls())
maf <- read.delim("data_mutations.txt")
maf$VAF <- maf$t_alt_count / (maf$t_alt_count + maf$t_ref_count)
maf$depth <- maf$t_alt_count + maf$t_ref_count
nsyv <- c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site')
maf <- subset(maf, Variant_Classification %in% nsyv)
genmut <- as.matrix( table(maf$Hugo_Symbol, maf$Tumor_Sample_Barcode) )
write.csv(genmut, file='MSK2017_mCRC_MSS_Gene_Level.csv',quote=F)
options(stringsAsFactors=F)
library(data.table)
setwd('<PROJECT_ROOT>')
rm(list=ls())
genmut <- read.csv('MSK2017_mCRC_MSS_Gene_Level.csv')
colnames(genmut)[1] <- 'symbol'
max(genmut[,2:ncol(genmut)])
for(i in 2:ncol(genmut)){
	genmut[,i] <- ifelse(genmut[,i]==0,0,1)
	rm(i)
}
max(genmut[,2:ncol(genmut)])
filename <- c('MSK2017_Early_stage_primaries',
			  'MSK2017_mCRC_Metastasis',
			  'MSK2017_mCRC_primaries',
			  'MSK2017_mCRC.PERITONEUM.OMENTUM.ABDOMEN_primaries',
			  'MSK2017_mCRC.PERITONEUM.OMENTUM.ABDOMEN_Metastasis')
for(grp in filename){
    sam <- read.delim(paste(grp,"_sampleID.txt",sep=''))
	sam$SampleID <- gsub('-','.',sam$SampleID)
	sda <- genmut[,intersect(sam$SampleID , colnames(genmut))]
	rownames(sda) <- genmut$symbol
	re <- data.frame()
    for(n in 1:nrow(sda)){
	sre <-  data.frame(gene = rownames(sda)[n],
		  Early_Stage_Primary_T = sum(sda[n,]),
		  Early_Stage_Primary_F = ncol(sda)-sum(sda[n,])
		  )
	re <- rbind(re,sre)
	rm(sre, n)
     }
	write.csv(re, file=paste(grp,"_Gene_Level_No.Mutation.csv",sep=''),row.names=F,quote=F)
	rm(sam,sda,grp,re)
}
