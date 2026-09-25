R
setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
mut <- read.delim('../results_mutDriver_paired/NJMU_mCRC_Nodule.maf')
maf <- read.delim('mCRC_candidate_driGene_mutation.maf')
maf <- subset(maf, Hugo_Symbol == 'SMARCB1' & paste(mutation_id.x,Normal,Class_a) %in% paste(mut$mutation_id,mut$Normal,mut$sample_id))
maf$vaf <- maf$t_alt_count/(maf$t_alt_count + maf$t_ref_count)
dat <- subset(maf, t_alt_count>=3 & t_alt_count/(t_alt_count + t_ref_count)>=0.01 & Class_a %in% c('Nodule-1','Nodule-2','Nodule-3','Nodule-4','Nodule-5','Nodule-6'))
maf <- subset(maf,!duplicated(paste(mutation_id.y, Normal)))
maf$ano <- NA
for(i in 1:nrow(maf)){
	ano <- strsplit(maf$Other_Transcripts[i],'\\|')[[1]]
	maf$ano[i] <- ano[grep('ENST00000263121',ano)]
	rm(ano,i)
}
maf$Annotation_Transcript <- sapply(strsplit(maf$ano,'_'),'[',2)
maf[,c('ano','Protein_Change','Variant_Classification')]
maf$Protein_Change <- c('p.E31V','p.R377H')
dat <- maf[,c('Normal','Hugo_Symbol','Annotation_Transcript','Chromosome','Start_Position','Tumor_Seq_Allele1','Tumor_Seq_Allele2','Protein_Change','Variant_Classification','Class_a','t_alt_count','t_ref_count')]
dat$Type <- 'NJMU mCRC Metastasis'
dat$aachange <- ifelse(dat$aachange=='','.',dat$aachange)
write.table(dat,file='SMARCB1/NJMU_mCRC_SMARCB1_forLolliplot.txt',row.names=F,quote=F,sep='\t')
options(stringsAsFactors=F)
library(data.table)
setwd('<PROJECT_ROOT>')
rm(list=ls())
ep <- read.delim('../crc_msk_2017/MSK2017_Early_stage_primaries_sampleID.txt')
mp <- read.delim('../crc_msk_2017/MSK2017_mCRC.PERITONEUM.OMENTUM.ABDOMEN_primaries_sampleID.txt')
mm <- read.delim('../crc_msk_2017/MSK2017_mCRC.PERITONEUM.OMENTUM.ABDOMEN_Metastasis_sampleID.txt')
ep$group <- 'MSK Early stage primaries'
mp$group <- 'MSK Abdomen mCRC primaries'
mm$group <- 'MSK Abdomen mCRC metastasis'
sam <- rbind(ep[,c('SampleID','group')], mp[,c('SampleID','group')], mm[,c('SampleID','group')])
sam <- subset(sam, !duplicated(SampleID))
maf <- read.delim("../crc_msk_2017<PROJECT_ROOT>_mutations.txt")
maf$VAF <- maf$t_alt_count / (maf$t_alt_count + maf$t_ref_count)
maf$depth <- maf$t_alt_count + maf$t_ref_count
nsyv <- c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site')
maf <- subset(maf, Variant_Classification %in% nsyv)
maf <- subset(maf, Hugo_Symbol == 'SMARCB1')
mut <- merge(sam, maf, by.x='SampleID', by.y='Tumor_Sample_Barcode')
dat <- mut[,c('SampleID','Hugo_Symbol','Transcript_ID','Chromosome','Start_Position','Tumor_Seq_Allele1','Tumor_Seq_Allele2','HGVSp_Short','Variant_Classification','group','t_alt_count','t_ref_count')]
dat$aachange <- ifelse(dat$aachange=='','.',dat$aachange)
write.table(dat,file='SMARCB1/MSK_SMARCB1_forLolliplot.txt',row.names=F,quote=F,sep='\t')
setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
ear <- read.delim('../coadread_tcga_pan_can_atlas_2018/TCGA_COADREAD_StageI-III_MSS_sampleID.txt')
ear$group <- 'TCGA Early Stage'
sam <- ear[,c('PATIENT_ID','group')]
maf <- fread('../coadread_tcga_pan_can_atlas_2018<PROJECT_ROOT>_mutations.txt',data.table=F)
nsyv <- c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site')
maf <- subset(maf, Variant_Classification %in% nsyv)
maf$VAF <- maf$t_alt_count / (maf$t_alt_count + maf$t_ref_count)
maf$depth <- maf$t_alt_count + maf$t_ref_count
maf$PATIENT_ID <- substr(maf$Tumor_Sample_Barcode,1,12)
maf <- subset(maf,  Hugo_Symbol == 'SMARCB1')
maf <- merge(sam, maf, by='PATIENT_ID')
dat <- maf[,c('PATIENT_ID','Hugo_Symbol','Transcript_ID','Chromosome','Start_Position','Tumor_Seq_Allele1','Tumor_Seq_Allele2','HGVSp_Short','Variant_Classification','group','t_alt_count','t_ref_count')]
dat$aachange <- ifelse(dat$aachange=='','.',dat$aachange)
write.table(dat,file='SMARCB1/TCGA_SMARCB1_forLolliplot.txt',row.names=F,quote=F,sep='\t')
options(stringsAsFactors=F)
setwd('<PROJECT_ROOT>')
library(data.table)
rm(list=ls())
maf <- read.delim('Zhongshan_candidate_genes.txt')
maf <- subset(maf, Hugo_Symbol == "SMARCB1")
dat <- maf[,c('Tumor_Sample_Barcode','Hugo_Symbol','Transcript_ID','Chromosome','Start_Position','Tumor_Seq_Allele1','Tumor_Seq_Allele2','HGVSp_Short','Variant_Classification','t_alt_count','t_ref_count')]
dat$refseq <- sapply(strsplit(dat$refseq,'\\.'),'[',1)
dat$aachange <- ifelse(dat$aachange=='','.',dat$aachange)
write.table(dat,file='SMARCB1/Zhongshan_SMARCB1_forLolliplot.txt',row.names=F,quote=F,sep='\t')
options(stringsAsFactors=F)
setwd('<PROJECT_ROOT>')
library(data.table)
rm(list=ls())
a <- read.delim('SMARCB1/NJMU_mCRC_SMARCB1_forLolliplot.txt')
b <- read.delim('SMARCB1/MSK_SMARCB1_forLolliplot.txt')
d <- read.delim('SMARCB1/TCGA_SMARCB1_forLolliplot.txt')
e <- read.delim('SMARCB1/Zhongshan_SMARCB1_forLolliplot.txt')
e$Type <- 'ZS Early Stage'
dat <- rbind(a,b,d,e)
dat$Note <- ifelse(dat$Type %in% c('MSK Abdomen mCRC metastasis','MSK Abdomen mCRC primaries','NJMU mCRC Metastasis'), 'mCRC', ifelse(dat$Type %in% c('MSK Early stage primaries','TCGA Early Stage','ZS Early Stage'), 'Early stage primary',NA ))
dat$refseq <- sapply(strsplit(dat$refseq,'\\.'),'[',1)
dat <- subset(dat,!duplicated(paste(Normal,chromosome,start,REF,ALT)))
dat$refseq <- sapply(strsplit(dat$refseq,'\\.'),'[',1)
write.table(dat,file='SMARCB1/SMARCB1_forLolliplot.txt',row.names=F,quote=F,sep='\t')
by(dat$Normal,dat$Type,function(x){length(unique(x))})
