cd <PROJECT_ROOT>
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
data <- data.frame()
for(n in unique(tre$Normal)){
  cln <- read.delim(paste(work_dir,'Pyclone/result/Input_all_VAF005_read4/',n, "/tables/loci.tsv",sep=''),comment.char = "#",sep='\t')
  clu_use <- read.delim(paste(n,'_cluster_use.txt',sep=''),h=F)
  cln <- subset(cln, cluster_id %in% clu_use$V1)
  cln$Class <- sapply(strsplit(cln$sample_id,'-'),'[',1)
  cln <- subset(cln, Class == 'Nodule')
  sda <- subset(cln, variant_allele_frequency > 0)
  sda$Normal <- n
  data <- rbind(data,sda)
  rm(sda, cln, n,clu_use)
}
write.table(data,file='NJMU_mCRC_Nodule.maf',row.names=F,quote=F,sep='\t')
system('cp NJMU_mCRC_Nodule.maf <PROJECT_ROOT>')
setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
maf <- fread('NJMU_mCRC_Nodule.maf',data.table=F)
maf$gene <- sapply(strsplit(maf$mutation_id,':'),'[',1)
maf$Variant_Classification <- sapply(strsplit(maf$mutation_id,':'),'[',2)
nsyv <- c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site')
maf <- subset(maf, Variant_Classification %in% nsyv)
genmut <- as.matrix( table(maf$gene, maf$Normal) )
write.csv(genmut, file='mCRC_MSS_Nodule_Gene_Level.csv',quote=F)
options(stringsAsFactors=F)
library(data.table)
setwd('<PROJECT_ROOT>')
rm(list=ls())
genmut <- read.csv('mCRC_MSS_Nodule_Gene_Level.csv')
colnames(genmut)[1] <- 'symbol'
max(genmut[,2:ncol(genmut)])
for(i in 2:ncol(genmut)){
	genmut[,i] <- ifelse(genmut[,i]==0,0,1)
	rm(i)
}
max(genmut[,2:ncol(genmut)])
sda <- genmut[,2:ncol(genmut)]
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
write.csv(re, file="mCRC_Gene_Level_No.Mutation.csv",row.names=F,quote=F)
