setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
info <- read.delim('data_clinical_patient.txt')
info <- subset(info, ! SUBTYPE %in% c('COAD_MSI','COAD_POLE','READ_MSI','READ_POLE'))
adv <- info
info <- subset(info, AJCC_PATHOLOGIC_TUMOR_STAGE %in% c('STAGE I','STAGE IA','STAGE II','STAGE IIA','STAGE IIB','STAGE IIC','STAGE III','STAGE IIIA','STAGE IIIB','STAGE IIIC'))
write.table(info, file='TCGA_COADREAD_StageI-III_MSS_sampleID.txt',row.names=F,quote=F,sep='\t')
setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
sam <- read.delim("TCGA_COADREAD_StageI-III_MSS_sampleID.txt")
sam$SampleID <- sam$PATIENT_ID
maf <- fread('data_mutations.txt',data.table=F)
maf$sampleID <- substr(maf$Tumor_Sample_Barcode,1,12)
nsyv <- c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site')
maf <- subset(maf, Variant_Classification %in% nsyv)
maf$VAF <- maf$t_alt_count / (maf$t_alt_count + maf$t_ref_count)
maf$depth <- maf$t_alt_count + maf$t_ref_count
maf$sampleID <- substr(maf$Tumor_Sample_Barcode,1,12)
genmut <- as.matrix( table(maf$Hugo_Symbol, maf$sampleID) )
write.csv(genmut, file='TCGA2018_COADREAD_EarlyStage_MSS_Gene_Level.csv',quote=F)
options(stringsAsFactors=F)
library(data.table)
setwd('<PROJECT_ROOT>')
rm(list=ls())
genmut <- read.csv('TCGA2018_COADREAD_EarlyStage_MSS_Gene_Level.csv')
colnames(genmut)[1] <- 'symbol'
max(genmut[,2:ncol(genmut)])
for(i in 2:ncol(genmut)){
	genmut[,i] <- ifelse(genmut[,i]==0,0,1)
	rm(i)
}
max(genmut[,2:ncol(genmut)])
filename <- c('TCGA_COADREAD_StageI-III_MSS')
for(grp in filename){
    sam <- read.delim(paste(grp,"_sampleID.txt",sep=''))
	sam$SampleID <- gsub('-','.',sam$PATIENT_ID)
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
