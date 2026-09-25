mkdir <PROJECT_ROOT>
cd <PROJECT_ROOT>
cp <PROJECT_ROOT> <PROJECT_ROOT>
cd <PROJECT_ROOT>
ls
R
setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
work_dir <- "<PROJECT_ROOT>"
info <- fread(paste(work_dir,'config/tumor_normal.list',sep=''),data.table=F)
info <- subset(info, Normal != 'CASE_ID')
for(i in 1:nrow(info)){
  maf <- read.csv(paste(work_dir,'maf/',info$Tumor[i],"_",info$Normal[i], "_GGA_Filter_funcotated.maf",sep=''),comment.char = "#",sep='\t')
  maf <- maf[,c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Strand','Variant_Classification','Variant_Type','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2','t_alt_count','t_ref_count','n_alt_count','n_ref_count')]
  write.table(maf, file=paste(info$Tumor[i],"_",info$Normal[i], "_GGA_Filter_use.txt",sep=''), row.names=F,quote=F,sep='\t')
  rm(maf)
}
tre <- read.csv('../Clonevol_sample_used.csv')
tre$Class <- sapply(strsplit(tre$Class_a,'-'),'[',1)
for(n in unique(tre$Normal)){
  cln <- read.delim(paste(work_dir,'Pyclone/result/Input_all_VAF005_read4/',n, "/tables/loci.tsv",sep=''),comment.char = "#",sep='\t')
  clu_use <- read.delim(paste(n,'_cluster_use.txt',sep=''),h=F)
  cln <- subset(cln, cluster_id %in% clu_use$V1)
  stre <- subset(tre, Normal == n & Class=='Primary')
  sinf <- subset(info, Normal == n & Class=='Primary' & Class_a %in% stre$Class_a)
  rm(stre)
  if(nrow(sinf)==1){
    scln <- cln[cln$sample_id==sinf$Class_a,]
	scln$clonality <- ifelse(scln$cellular_prevalence>0.5,'clonal',ifelse(scln$cellular_prevalence>0.1,'subclonal','absent'))
	write.table(scln, file=paste(n,'_primary_clonality_0.5_0.1.tsv',sep=''), row.names=F, quote=F, sep='\t')
	rm(scln)
    }
  if(nrow(sinf)==2){
    scln_1 <- cln[cln$sample_id==sinf$Class_a[1],]
	maf_1 <- read.delim(paste(sinf[1,'Tumor'],"_",sinf[1,'Normal'], "_GGA_Filter_use.txt",sep=''))
	maf_1$mutation_id <- paste(maf_1$Hugo_Symbol,maf_1$Variant_Classification,maf_1$Chromosome,maf_1$Start_Position,maf_1$Tumor_Seq_Allele1,maf_1$Tumor_Seq_Allele2,sep=':')
	maf_1$depth <- maf_1$t_alt_count + maf_1$t_ref_count
	pri_1 <- merge(scln_1, maf_1[,c('mutation_id','depth')], by='mutation_id')
    scln_2 <- cln[cln$sample_id==sinf$Class_a[2],]
	maf_2 <- read.delim(paste(sinf[2,'Tumor'],"_",sinf[2,'Normal'], "_GGA_Filter_use.txt",sep=''))
	maf_2$mutation_id <- paste(maf_2$Hugo_Symbol,maf_2$Variant_Classification,maf_2$Chromosome,maf_2$Start_Position,maf_2$Tumor_Seq_Allele1,maf_2$Tumor_Seq_Allele2,sep=':')
	maf_2$depth <- maf_2$t_alt_count + maf_2$t_ref_count
	pri_2 <- merge(scln_2, maf_2[,c('mutation_id','depth')], by='mutation_id')
	scln <- merge(pri_1,pri_2,by='mutation_id',all=T)
	scln$CCF <- (scln$cellular_prevalence.p1*scln$depth.p1 + scln$cellular_prevalence.p2*scln$depth.p2)/(scln$depth.p1 + scln$depth.p2)
	scln$clonality <- ifelse(scln$CCF>0.5,'clonal',ifelse(scln$CCF>0.1,'subclonal','absent'))
	write.table(scln, file=paste(n,'_primary_clonality_0.5_0.1.tsv',sep=''), row.names=F, quote=F, sep='\t')
	rm(scln, scln_1, maf_1, pri_1, scln_2, maf_2, pri_2)
    }
  if(nrow(sinf)==3){
    scln_1 <- cln[cln$sample_id==sinf$Class_a[1],]
	maf_1 <- read.delim(paste(sinf[1,'Tumor'],"_",sinf[1,'Normal'], "_GGA_Filter_use.txt",sep=''))
	maf_1$mutation_id <- paste(maf_1$Hugo_Symbol,maf_1$Variant_Classification,maf_1$Chromosome,maf_1$Start_Position,maf_1$Tumor_Seq_Allele1,maf_1$Tumor_Seq_Allele2,sep=':')
	maf_1$depth <- maf_1$t_alt_count + maf_1$t_ref_count
	pri_1 <- merge(scln_1, maf_1[,c('mutation_id','depth')], by='mutation_id')
    scln_2 <- cln[cln$sample_id==sinf$Class_a[2],]
	maf_2 <- read.delim(paste(sinf[2,'Tumor'],"_",sinf[2,'Normal'], "_GGA_Filter_use.txt",sep=''))
	maf_2$mutation_id <- paste(maf_2$Hugo_Symbol,maf_2$Variant_Classification,maf_2$Chromosome,maf_2$Start_Position,maf_2$Tumor_Seq_Allele1,maf_2$Tumor_Seq_Allele2,sep=':')
	maf_2$depth <- maf_2$t_alt_count + maf_2$t_ref_count
	pri_2 <- merge(scln_2, maf_2[,c('mutation_id','depth')], by='mutation_id')
    scln_3 <- cln[cln$sample_id==sinf$Class_a[3],]
	maf_3 <- read.delim(paste(sinf[3,'Tumor'],"_",sinf[3,'Normal'], "_GGA_Filter_use.txt",sep=''))
	maf_3$mutation_id <- paste(maf_3$Hugo_Symbol,maf_3$Variant_Classification,maf_3$Chromosome,maf_3$Start_Position,maf_3$Tumor_Seq_Allele1,maf_3$Tumor_Seq_Allele2,sep=':')
	maf_3$depth <- maf_3$t_alt_count + maf_3$t_ref_count
	pri_3 <- merge(scln_3, maf_3[,c('mutation_id','depth')], by='mutation_id')
	scln <- merge(pri_1,pri_2,by='mutation_id',all=T)
	scln <- merge(scln,pri_3,by='mutation_id',all=T)
	scln$CCF <- (scln$cellular_prevalence.p1*scln$depth.p1 + scln$cellular_prevalence.p2*scln$depth.p2 + scln$cellular_prevalence.p3*scln$depth.p3)/(scln$depth.p1 + scln$depth.p2 + scln$depth.p3)
	scln$clonality <- ifelse(scln$CCF>0.5,'clonal',ifelse(scln$CCF>0.1,'subclonal','absent'))
	write.table(scln, file=paste(n,'_primary_clonality_0.5_0.1.tsv',sep=''), row.names=F, quote=F, sep='\t')
	rm(scln, scln_1, maf_1, pri_1, scln_2, maf_2, pri_2, scln_3, maf_3, pri_3)
    }
  rm(sinf, cln, n)
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
  pri <- read.delim(paste(n,'_primary_clonality_0.5_0.1.tsv',sep=''))
  if(!'CCF'%in%colnames(pri)){pri$CCF <- pri$cellular_prevalence}
  maf <- read.delim(paste(work_dir,'Pyclone/result/Input_all_VAF005_read4/',n, "/tables/loci.tsv",sep=''),comment.char = "#",sep='\t')
  clu_use <- read.delim(paste(n,'_cluster_use.txt',sep=''),h=F)
  maf <- merge(maf,clu_use,by='cluster_id')
  stre <- subset(tre, Normal == n & Class=='Nodule')
  maf <- subset(maf, sample_id %in% stre$Class_a)
  maf$clonality <- ifelse(maf$cellular_prevalence>0.5,'clonal',ifelse(maf$cellular_prevalence>0.1,'subclonal','absent'))
  colnames(maf)[ncol(maf)] <- 'clonality.met'
  sda <- merge(pri[,c('mutation_id','clonality','CCF')], maf, by='mutation_id')
  sda$Normal = n
  data <- rbind(data, sda)
  rm(pri, maf, sda, stre, n, clu_use)
}
data$clone <- paste(data$clonality, data$clonality.met, sep='-')
data$type <- ifelse( data$clone %in% c('absent-clonal','absent-subclonal','subclonal-clonal'), "metastasis favored", ifelse( data$clone %in% c('clonal-absent','subclonal-absent'), "primary favored", ifelse(data$clone %in% c('clonal-clonal','subclonal-subclonal','clonal-subclonal'), "maintained", NA) ))
sum(table(data$type))
write.csv(data, file='NJMU_mCRC_mutation_clonality_primary_metastases_vaf005_read4_0.5_0.1.csv',row.names=F,quote=F)
R
setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
data <- read.csv("NJMU_mCRC_mutation_clonality_primary_metastases_vaf005_read4_0.5_0.1.csv")
dat1 <- subset(data, ! Normal %in% c('CASE_ID','CASE_ID'))
filelist <- data.frame(filename=c(
'CASE_ID_Nodule-1.tsv',
'CASE_ID_Nodule-2.tsv',
'CASE_ID_Nodule-1.tsv',
'CASE_ID_Nodule-2.tsv',
'CASE_ID_Nodule-3.tsv'
),
Class = c(
'Nodule-1',
'Nodule-2',
'Nodule-1',
'Nodule-2',
'Nodule-3'
))
filelist
dat2 <- data.frame()
for(i in 1:nrow(filelist)){
    snorm <- subset(data, Normal %in% strsplit(filelist$filename[i],'_')[[1]][1])
	snorm <- subset(snorm, !duplicated(mutation_id))
	snorm <- snorm[,c('mutation_id','clonality.primary','CCF','cluster_id')]
	sda <- read.delim(paste('../../9_time_dissemination_NJMU/',filelist$filename[i],sep=''))
	sda <- sda[,c('mutation_id','CCF')]
	colnames(sda)[2] <- 'cellular_prevalence'
	sda$sample_id <- filelist$Class[i]
	sda$clonality.metastases <- ifelse(sda$cellular_prevalence>0.5,'clonal',ifelse(sda$cellular_prevalence>0.1,'subclonal','absent'))
	sda$Normal <- strsplit(filelist$filename[i],'_')[[1]][1]
	sda1 <- merge(snorm,sda,by='mutation_id')
	dat2 <- rbind(dat2, sda1)
	rm(snorm, sda, sda1)
}
dat2$clone <- paste(dat2$clonality.primary, dat2$clonality.metastases, sep='-')
dat2$type <- ifelse( dat2$clone %in% c('absent-clonal','absent-subclonal','subclonal-clonal'), "metastasis favored", ifelse( dat2$clone %in% c('clonal-absent','subclonal-absent'), "primary favored", ifelse(dat2$clone %in% c('clonal-clonal','subclonal-subclonal','clonal-subclonal'), "maintained", NA) ))
dat <- rbind(dat1[,colnames(dat1)%in%colnames(dat2)],dat2)
write.csv(dat, file='NJMU_mCRC_mutation_clonality_primary_metastases_vaf005_read4_0.5_0.1_corrected.csv',row.names=F,quote=F)
setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
data <- read.csv("NJMU_mCRC_mutation_clonality_primary_metastases_vaf005_read4_0.5_0.1_corrected.csv")
data <- subset(data, clone != 'absent-absent')
data$Variant_Classification <- sapply(strsplit(data$mutation_id,':'),'[',2)
data$gene <- sapply(strsplit(data$mutation_id,':'),'[',1)
nsyv <- c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site')
drimut <- subset(data, Variant_Classification %in% nsyv)
genmut <- subset(drimut, clonality.metastases != 'absent')
genmut <- subset(genmut , cellular_prevalence>0.5 | (cellular_prevalence-CCF>0.1))
count_sample <- function(x){ return ( length(unique(x)) ) }
re <- as.matrix (by(genmut$Normal, genmut$gene, count_sample) )
re <- as.data.frame(re)
colnames(re) <- 'n_sample'
re$gene <- rownames(re)
dri <- read.csv('<PROJECT_ROOT>')
dri <- subset(dri, database %in% c('IntOGen-COREAD','TCGA_Nature.2012'))
cosmic <- read.csv('<PROJECT_ROOT>',h=F)
nodri <- read.csv('<PROJECT_ROOT>')
cosmic <- subset(cosmic, !V1 %in% nodri$GeneSymbol)
tcga <- read.csv('<PROJECT_ROOT>')
msk <- read.table('<PROJECT_ROOT>')
driver <- unique(c(dri$SYMBOL, cosmic$V1, tcga$Gene, msk$V1))
drim <- subset(re, gene %in% driver)
cand <- subset(re, n_sample >= 2)
gen_include <- rbind(drim, cand)
gen_include <- subset(gen_include,!duplicated(gene))
write.csv(gen_include, file='NJMU_mCRC_Included_genes_0.5_0.1_corrected.csv',row.names=F,quote=F)
system('cp NJMU_mCRC_Included_genes_0.5_0.1_corrected.csv <PROJECT_ROOT>')
dri <- read.csv('../Cancer_Mut_Driver_Genes.csv')
drimat <- as.matrix(table(drimut$gene, drimut$type))
write.table(drimat,file='NJMU_All_Gene_Mutation_Type_vaf005_0.5_0.1_corrected.txt',row.names=T,quote=F,sep='\t')
drimat <- read.delim('NJMU_All_Gene_Mutation_Type_vaf005_0.5_0.1_corrected.txt')
drimat$gene <- rownames(drimat)
drimat <- subset(drimat, gene %in% c(drim$gene,cand$gene))
simut <- subset(data, (!gene %in% dri$SYMBOL) & !Variant_Classification %in% nsyv)
prob <- as.numeric(table(simut$type)/sum(table(simut$type)))
library(EMT)
multinomial.Chisq <- function(x){
  test <- multinomial.test( x, c(0.3191863, 0.3103968, 0.3704169), useChisq = TRUE )
  return(test$p.value)
  }
drimat$multinomial.chisq <- apply(drimat[,1:3],1, multinomial.Chisq)
write.csv(drimat,file='NJMU_All.Metastasis.Driver.Gene_Paired_Analysis_Results_vaf005_0.5_0.1_corrected.csv',row.names=F,quote=F)
rm(list=ls())
dat <- read.csv('NJMU_All.Metastasis.Driver.Gene_Paired_Analysis_Results_vaf005_0.5_0.1_corrected.csv')
dat$no.mut <- dat$maintained + dat$metastasis.favored + dat$primary.favored
dat <- subset(dat, no.mut>=5)
dat$metastasis.prop <- dat$metastasis.favored / dat$no.mut
dat$maintained.prop <- dat$maintained / dat$no.mut
dat$primary.prop <- dat$primary.favored / dat$no.mut
dat$logFC.metastasis <- log2(dat$metastasis.prop /  0.3103968 )
dat$logFC.maintained <- log2(dat$maintained.prop / 0.3191863)
dat$logFC.primary <- log2(dat$primary.prop / 0.3704169)
dat$multinomial.chisq.FDR <- p.adjust(dat$multinomial.chisq, method='BH')
subset(dat, gene %in% c('TP53','KRAS','APC','SMAD4'))
subset(dat, gene %in% c('RNF43','SMARCB1','BRAF','PTEN'))
sig <- subset(dat, multinomial.chisq<0.05)
write.csv(dat,file='NJMU_mCRC_Metastasis_Driver_Gene_Paired_Analysis_Results_vaf005_0.5_0.1.csv',row.names=F,quote=F)
write.csv(sig,file='NJMU_mCRC_Metastasis_Driver_Gene_Paired_Analysis_P0.05_vaf005_0.5_0.1.csv',row.names=F,quote=F)
subset(sig, gene %in% c('TP53','KRAS','APC','SMAD4'))
system('cp NJMU_mCRC_Metastasis_Driver_Gene_Paired_Analysis_P0.05_vaf005_0.5_0.1.csv <PROJECT_ROOT>')
system('cp NJMU_mCRC_Metastasis_Driver_Gene_Paired_Analysis_Results_vaf005_0.5_0.1.csv <PROJECT_ROOT>')
