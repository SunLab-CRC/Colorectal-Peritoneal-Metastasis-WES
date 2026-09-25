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
for(n in unique(tre$Normal)){
  cln <- read.delim(paste(work_dir,'Pyclone/result/Input_all_VAF005_read4/',n, "/tables/loci.tsv",sep=''),comment.char = "#",sep='\t')
  clu_use <- read.delim(paste('../6_Met_driver_genes/clone_use_vaf005_0801/',n,'_cluster_use.txt',sep=''),h=F)
  cln <- subset(cln, cluster_id %in% clu_use$V1)
  cln$Class <- sapply(strsplit(cln$sample_id,'-'),'[',1)
  sda <- subset(cln, variant_allele_frequency > 0)
  sda$Normal <- n
  data <- rbind(data,sda)
  rm(sda, cln, n,clu_use)
}
write.table(data,file='NJMU_mCRC_ALL_Sample_mutation.maf',row.names=F,quote=F,sep='\t')
system('cp NJMU_mCRC_ALL_Sample_mutation.maf <PROJECT_ROOT>')
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
info <- subset(info, paste(Normal, Class_a) %in% paste(tre$Normal, tre$Class_a) )
nsyv <- c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site')
data <- data.frame()
for(normal in unique(info$Normal)){
  sinf <- subset(info, Normal == normal)
  sinf$Class_a <- gsub('-','_',sinf$Class_a)
  for(i in 1:nrow(sinf)){
    gga <- read.csv(paste('<PROJECT_ROOT>',sinf$Tumor[i],'_',sinf$Normal[i],'_GGA_Filter_funcotated.maf',sep=''), comment.char = "#", sep='\t')
    gga <- subset(gga,Variant_Classification %in% nsyv)
    gga$mutID <- paste(gga$Hugo_Symbol, gga$Variant_Classification, gga$Chromosome, gga$Start_Position, gga$Reference_Allele, gga$Tumor_Seq_Allele2, sep=':')
	dat <- gga[,c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2','Genome_Change','Annotation_Transcript','Transcript_Strand','Transcript_Exon','Transcript_Position','cDNA_Change','Codon_Change','Protein_Change','Refseq_mRNA_Id','Refseq_prot_Id','t_alt_count','t_ref_count','HGNC_HGNC_ID','HGNC_Chromosome','HGNC_Ensembl_Gene_ID','HGNC_RefSeq_IDs','HGNC_RefSeq.supplied_by_NCBI.','mutID','Other_Transcripts')]
	dat$Normal <- sinf$Normal[i]
	dat$Tumor <- sinf$Tumor[i]
	dat$Class_a <- sinf$Class_a[i]
	dat$Class <- sinf$Class[i]
    data <- rbind(data, dat)
    rm(gga, dat, i)
    }
rm(sinf)
rm(normal)
}
write.table(data,file='NJMU_mCRC_ALL_Sample_mutation_annotation.maf',row.names=F,quote=F,sep='\t')
system('cp NJMU_mCRC_ALL_Sample_mutation_annotation.maf <PROJECT_ROOT>')
