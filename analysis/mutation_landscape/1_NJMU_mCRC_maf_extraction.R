mkdir <PROJECT_ROOT>
cd <PROJECT_ROOT>
R
setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
work_dir <- "<PROJECT_ROOT>"
info <- fread(paste(work_dir,'config/tumor_normal.list',sep=''),data.table=F)
info <- subset(info, Normal != 'CASE_ID')
tre <- read.csv('../6_Met_driver_genes/Clonevol_sample_used.csv')
tre$Class <- sapply(strsplit(tre$Class_a,'-'),'[',1)
data <- data.frame()
for(n in unique(info$Normal)){
  sinf <- subset(info,Normal==n)
  stre <- subset(tre,Normal==n)
  sinf <- subset(sinf, Class_a %in% stre$Class_a)
  rm(stre)
  maf <- data.frame()
  for(i in 1:nrow(sinf)){
    smaf <- read.csv(paste(work_dir,'maf/',sinf$Tumor[i],"_",sinf$Normal[i], "_GGA_Filter_funcotated.maf",sep=''),comment.char = "#",sep='\t')
	smaf$Class_a <- sinf$Class_a[i]
	maf <- rbind(maf, smaf)
	rm(smaf, i)
  }
  maf$mutation_id <- paste(maf$Hugo_Symbol,maf$Variant_Classification,maf$Chromosome,maf$Start_Position,maf$Tumor_Seq_Allele1,maf$Tumor_Seq_Allele2,sep=':')
  cln <- read.delim(paste(work_dir,'Pyclone/result/Input_all_VAF005_read4/',n, "/tables/loci.tsv",sep=''),comment.char = "#",sep='\t')
   clu_use <- read.delim(paste('../6_Met_driver_genes/clone_use_vaf005_0801/',n,'_cluster_use.txt',sep=''),h=F)
   cln <- merge(cln, clu_use,by='cluster_id')
   cln$id <- paste(cln$mutation_id, cln$sample_id, sep='_')
   maf$id <- paste(maf$mutation_id, maf$Class_a, sep='_')
   sda <- merge(cln,maf,by='id')
   sda$Normal = n
   data <- rbind(data,sda)
  rm(sinf, maf, clu_use, cln, sda, n)
}
write.table(data,file='NJMU_mCRC_ALL.Gene_mutation.maf',row.names=F,quote=F,sep='\t')
system('sz NJMU_mCRC_ALL.Gene_mutation.maf')
