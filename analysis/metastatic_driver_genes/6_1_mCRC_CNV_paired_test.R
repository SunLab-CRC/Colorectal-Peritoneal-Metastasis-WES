cd <PROJECT_ROOT>
cp <PROJECT_ROOT> ./
<PROJECT_ROOT> intersect \
-a SCNA_driver_crc/gencode.v19.annotation.bed \
-b cytoBand.txt \
-wa -wb \
> gencode.v19.annotation_cytoband.bed
wc -l gencode.v19.annotation_cytoband.bed
head <PROJECT_ROOT>
R
setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
rm(list=ls())
segs <- read.delim('<PROJECT_ROOT>')
segs$chr <- paste('chr', segs$Chromosome, sep='')
by(segs$Corrected_Copy_Number, segs$Corrected_Call, summary)
dat <- segs[,c('chr','Start','End','Sample','Corrected_Call')]
dat <- subset(dat, Corrected_Call %in% c('AMP','HLAMP','HETD','HOMD'))
write.table(dat, file='mCRC_Titan_all_seg.final.bed',row.names=F,col.names=F,quote=F,sep='\t')
cd <PROJECT_ROOT>
wc -l mCRC_Titan_all_seg.final.bed
<PROJECT_ROOT> intersect -a ../SCNA_driver_crc/gencode.v19.annotation.bed -b mCRC_Titan_all_seg.final.bed -wa -wb -f 0.5 > mCRC_SCNA_gene_sample.bed
wc -l mCRC_SCNA_gene_sample.bed
cp mCRC_SCNA_gene_sample.bed <PROJECT_ROOT>
R
setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
work_dir <- "<PROJECT_ROOT>"
info <- fread(paste(work_dir,'config/tumor_normal.list',sep=''),data.table=F)
info <- subset(info, Normal != 'CASE_ID')
tre <- read.csv('../Clonevol_sample_used.csv')
tre$Class <- sapply(strsplit(tre$Class_a,'-'),'[',1)
for(n in unique(tre$Normal)){
  cln <- read.delim('mCRC_SCNA_gene_sample.bed',h=F)
  cln <- subset(cln, V7 == ' protein_coding')
  cln$Normal <- sapply(strsplit(cln$V11,'_'),'[',2)
  cln$Tumor <- sapply(strsplit(cln$V11,'_'),'[',1)
  stre <- subset(tre, Normal == n & Class=='Primary')
  sinf <- subset(info, Normal == n & Class=='Primary' & Class_a %in% stre$Class_a)
  rm(stre)
  scln <- subset(cln, Tumor %in% sinf$Tumor & Normal %in% sinf$Normal)
  scln$gene <- paste(scln$V5, scln$V6, sep='')
  pri <- as.matrix(table(scln$gene, scln$V12))
  write.table(pri, file=paste(n,'primary.txt',sep=''),quote=F,sep='\t')
  prid <- read.delim( paste(n,'primary.txt',sep='') )
  prid$amp <- apply(prid[,colnames(prid) %in% c('AMP','HLAMP'),drop=F],1,sum)
  prid$del <- apply(prid[,colnames(prid) %in% c('HETD','HOMD'),drop=F],1,sum)
  prid$gene <- rownames(prid)
  prid <- prid[,c('gene','amp','del')]
  rm(scln, pri, sinf)
  stre <- subset(tre, Normal == n & Class=='Nodule')
  sinf <- subset(info, Normal == n & Class=='Nodule' & Class_a %in% stre$Class_a)
  rm(stre)
  scln <- subset(cln, Tumor %in% sinf$Tumor & Normal %in% sinf$Normal)
  scln$gene <- paste(scln$V5, scln$V6, sep='')
  scln$sam_gene <- paste(scln$V11, scln$gene, sep=':')
  met <- as.matrix(table(scln$sam_gene, scln$V12))
  write.table(met, file=paste(n,'Nodule.txt',sep=''),quote=F,sep='\t')
  metd <- read.delim( paste(n,'Nodule.txt',sep='') )
  metd$amp <- apply(metd[,colnames(metd) %in% c('AMP','HLAMP'),drop=F],1,sum)
  metd$del <- apply(metd[,colnames(metd) %in% c('HETD','HOMD'),drop=F],1,sum)
  metd$sam_gene <- rownames(metd)
  metd$gene <- sapply(strsplit(metd$sam_gene,':'),'[',2)
  metd <- metd[,c('sam_gene','gene','amp','del')]
  sre <- merge(prid, metd, by='gene', all=T)
  a <- subset(sre, is.na(AMP.primary))
  a$AMP.primary <- 0
  a$DEL.primary <- 0
  b <- data.frame()
  for(k in unique(scln$V11)){
   sb <- subset(sre, is.na(AMP.Nodule))
   if(nrow(sb) == 0){next}
   sb$sam_gene <- paste(k, sb$gene, sep=':')
   sb$AMP.Nodule <- 0
   sb$DEL.Nodule <- 0
   b <- rbind(b,sb)
   rm(sb,k)
  }
  d <- subset(sre, !is.na(AMP.primary) & !is.na(AMP.Nodule))
  sred <- rbind(a,d,b)
  write.table(sred, file=paste(n,'_SCNA_gene_result.txt',sep=''),row.names=F, quote=F,sep='\t')
  rm(scln, met, sinf, sre, a, b, sred, n, prid, cln, d, metd)
    }
setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
work_dir <- "<PROJECT_ROOT>"
info <- fread(paste(work_dir,'config/tumor_normal.list',sep=''),data.table=F)
info <- subset(info, Normal != 'CASE_ID')
tre <- read.csv('../Clonevol_sample_used.csv')
tre$Class <- sapply(strsplit(tre$Class_a,'-'),'[',1)
data <- data.frame()
for(n in unique(info$Normal)){
  sda <- read.delim(paste(n,'_SCNA_gene_result.txt',sep=''))
  data <- rbind(data, sda)
  rm(sda, n)
}
data$samID <- sapply(strsplit(data$sam_gene,':'),'[',1)
data$AMP.primary <- ifelse(data$AMP.primary==0,0,1)
data$DEL.primary <- ifelse(data$DEL.primary==0,0,1)
data$AMP.Nodule <- ifelse(data$AMP.Nodule==0,0,1)
data$DEL.Nodule <- ifelse(data$DEL.Nodule==0,0,1)
dat <- subset(data, ! ((AMP.primary==1 & DEL.primary==1) | (AMP.Nodule==1 & DEL.Nodule==1) ) )
data <- dat
data$AMP <- paste(data$AMP.primary, data$AMP.Nodule, sep='-')
amp <- subset(data, AMP != "0-0")
amp$type <- ifelse( amp$AMP == "1-1", "maintained", ifelse( amp$AMP == "1-0", "primary favored", ifelse( amp$AMP == "0-1", "metastasis favored", NA )  ))
write.csv(amp, file='NJMU_mCRC_SCNA_AMP_primary_metastases.csv',row.names=F,quote=F)
data$DEL <- paste(data$DEL.primary, data$DEL.Nodule, sep='-')
DEL <- subset(data, DEL != "0-0")
DEL$type <- ifelse( DEL$DEL == "1-1", "maintained", ifelse( DEL$DEL == "1-0", "primary favored", ifelse( DEL$DEL == "0-1", "metastasis favored", NA )  ))
write.csv(DEL, file='NJMU_mCRC_SCNA_DEL_primary_metastases.csv',row.names=F,quote=F)
setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
data <- read.csv("NJMU_mCRC_SCNA_AMP_primary_metastases.csv")
data$symbol <- sapply(strsplit(data$gene,' '),'[',2)
data$ID <- sapply(strsplit(data$gene,' '),'[',1)
data$ID <- sapply(strsplit(data$gene,'\\.'),'[',1)
data$symbol <- toupper(gsub(' ','',data$symbol))
data$ID <- toupper(gsub(' ','',data$ID))
dri <- read.delim('../SCNA_driver_crc/SCNA_driver_CRC_GENCODE.V19.bed',h=F)
dri <- subset(dri, V8 == 'AMP')
dat <- subset(data, (symbol %in% dri$V6) | (ID %in% dri$V5) )
drimat <- as.matrix(table(dat$gene, dat$type))
write.table(drimat,file='NJMU_driGene_AMP_Type.txt',row.names=T,quote=F,sep='\t')
drimat <- read.delim('NJMU_driGene_AMP_Type.txt')
drimat$gene <- rownames(drimat)
simut <- subset(data, (! symbol %in% dri$V6) & (! ID %in% dri$V5))
library(EMT)
multinomial.Chisq <- function(x){
  test <- multinomial.test( x, c(0.3118664, 0.2741716, 0.4139620), useChisq = TRUE )
  return(test$p.value)
  }
drimat$multinomial.chisq <- apply(drimat[,1:3],1, multinomial.Chisq)
write.csv(drimat,file='NJMU_Amplification.Metastasis.Driver.Gene_Paired_Analysis_Results.csv',row.names=F,quote=F)
rm(list=ls())
dat <- read.csv('NJMU_Amplification.Metastasis.Driver.Gene_Paired_Analysis_Results.csv')
gen <- read.table('../SCNA_driver_crc/gencode.v19.annotation.bed')
dat$ID <- sapply(strsplit(dat$gene,' '),'[',1)
dat1 <- merge(dat,gen,by.x='ID',by.y='V5')
band <- read.table('<PROJECT_ROOT>')
band$cytoBand <- paste(gsub('chr','',band$V1),band$V11,sep='')
dat2 <- merge(dat1,band[,c('V5','cytoBand')],by.x='ID',by.y='V5',all.x=T)
write.csv(dat2,file='NJMU_Amplification.Metastasis.Driver.Gene_Paired_Analysis_Results_cytoband.csv',row.names=F,quote=F)
system('cp NJMU_Amplification.Metastasis.Driver.Gene_Paired_Analysis_Results_cytoband.csv <PROJECT_ROOT>')
setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
data <- read.csv("NJMU_mCRC_SCNA_DEL_primary_metastases.csv")
data$symbol <- sapply(strsplit(data$gene,' '),'[',2)
data$ID <- sapply(strsplit(data$gene,' '),'[',1)
data$ID <- sapply(strsplit(data$gene,'\\.'),'[',1)
data$symbol <- toupper(gsub(' ','',data$symbol))
data$ID <- toupper(gsub(' ','',data$ID))
dri <- read.delim('../SCNA_driver_crc/SCNA_driver_CRC_GENCODE.V19.bed',h=F)
dri <- subset(dri, V8 == 'DEL')
dat <- subset(data, (symbol %in% dri$V6) | (ID %in% dri$V5) )
drimat <- as.matrix(table(dat$gene, dat$type))
write.table(drimat,file='NJMU_driGene_DEL_Type.txt',row.names=T,quote=F,sep='\t')
drimat <- read.delim('NJMU_driGene_DEL_Type.txt')
drimat$gene <- rownames(drimat)
simut <- subset(data, (! symbol %in% dri$V6) & (! ID %in% dri$V5))
library(EMT)
multinomial.Chisq <- function(x){
  test <- multinomial.test( x, c(0.5556084, 0.2849778, 0.1594138), useChisq = TRUE )
  return(test$p.value)
  }
drimat$multinomial.chisq <- apply(drimat[,1:3],1, multinomial.Chisq)
write.csv(drimat,file='NJMU_Deletion.Metastasis.Driver.Gene_Paired_Analysis_Results.csv',row.names=F,quote=F)
rm(list=ls())
dat <- read.csv('NJMU_Deletion.Metastasis.Driver.Gene_Paired_Analysis_Results.csv')
gen <- read.table('../SCNA_driver_crc/gencode.v19.annotation.bed')
dat$ID <- sapply(strsplit(dat$gene,' '),'[',1)
dat1 <- merge(dat,gen,by.x='ID',by.y='V5')
band <- read.table('<PROJECT_ROOT>')
band$cytoBand <- paste(gsub('chr','',band$V1),band$V11,sep='')
dat2 <- merge(dat1,band[,c('V5','cytoBand')],by.x='ID',by.y='V5',all.x=T)
write.csv(dat2,file='NJMU_Deletion.Metastasis.Driver.Gene_Paired_Analysis_Results_cytoband.csv',row.names=F,quote=F)
system('cp NJMU_Deletion.Metastasis.Driver.Gene_Paired_Analysis_Results_cytoband.csv <PROJECT_ROOT>')
