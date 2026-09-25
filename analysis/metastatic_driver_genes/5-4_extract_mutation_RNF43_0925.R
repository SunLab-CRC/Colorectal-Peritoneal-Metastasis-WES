R
setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
mut <- read.delim('../results_mutDriver_paired/NJMU_mCRC_Nodule.maf')
maf <- read.delim('mCRC_candidate_driGene_mutation.maf')
maf <- subset(maf, Hugo_Symbol == 'RNF43' & paste(mutation_id.x,Normal,Class_a) %in% paste(mut$mutation_id,mut$Normal,mut$sample_id))
maf$vaf <- maf$t_alt_count/(maf$t_alt_count + maf$t_ref_count)
dat <- subset(maf, t_alt_count>=3 & t_alt_count/(t_alt_count + t_ref_count)>=0.01 & Class_a %in% c('Nodule-1','Nodule-2','Nodule-3','Nodule-4','Nodule-5','Nodule-6'))
maf <- subset(maf,!duplicated(paste(mutation_id.y, Normal)))
maf$ano <- NA
for(i in 1:nrow(maf)){
	ano <- strsplit(maf$Other_Transcripts[i],'\\|')[[1]]
	maf$ano[i] <- ano[grep('ENST00000407977',ano)]
	rm(ano,i)
}
maf$Annotation_Transcript <- sapply(strsplit(maf$ano,'_'),'[',2)
dat <- maf[,c('Normal','Hugo_Symbol','Annotation_Transcript','Chromosome','Start_Position','Tumor_Seq_Allele1','Tumor_Seq_Allele2','Protein_Change','Variant_Classification','Class_a','t_alt_count','t_ref_count')]
dat$Type <- 'NJMU mCRC Metastasis'
dat$aachange <- ifelse(dat$aachange=='','.',dat$aachange)
write.table(dat,file='lolliplot/NJMU_mCRC_RNF43_forLolliplot.txt',row.names=F,quote=F,sep='\t')
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
maf <- subset(maf, Hugo_Symbol == 'RNF43')
mut <- merge(sam, maf, by.x='SampleID', by.y='Tumor_Sample_Barcode')
dat <- mut[,c('SampleID','Hugo_Symbol','Transcript_ID','Chromosome','Start_Position','Tumor_Seq_Allele1','Tumor_Seq_Allele2','HGVSp_Short','Variant_Classification','group','t_alt_count','t_ref_count')]
dat$aachange <- ifelse(dat$aachange=='','.',dat$aachange)
write.table(dat,file='lolliplot/MSK_RNF43_forLolliplot.txt',row.names=F,quote=F,sep='\t')
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
maf <- subset(maf,  Hugo_Symbol == 'RNF43')
maf <- merge(sam, maf, by='PATIENT_ID')
dat <- maf[,c('PATIENT_ID','Hugo_Symbol','Transcript_ID','Chromosome','Start_Position','Tumor_Seq_Allele1','Tumor_Seq_Allele2','HGVSp_Short','Variant_Classification','group','t_alt_count','t_ref_count')]
dat$aachange <- ifelse(dat$aachange=='','.',dat$aachange)
write.table(dat,file='lolliplot/TCGA_RNF43_forLolliplot.txt',row.names=F,quote=F,sep='\t')
options(stringsAsFactors=F)
setwd('<PROJECT_ROOT>')
library(data.table)
rm(list=ls())
a <- read.delim('lolliplot/NJMU_mCRC_RNF43_forLolliplot.txt')
b <- read.delim('lolliplot/MSK_RNF43_forLolliplot.txt')
d <- read.delim('lolliplot/TCGA_RNF43_forLolliplot.txt')
dat <- rbind(a,b,d)
dat$Note <- ifelse(dat$Type %in% c('MSK Abdomen mCRC metastasis','MSK Abdomen mCRC primaries','NJMU mCRC Metastasis'), 'mCRC', ifelse(dat$Type %in% c('MSK Early stage primaries','TCGA Early Stage','ZS Early Stage'), 'Early stage primary',NA ))
dat$refseq <- sapply(strsplit(dat$refseq,'\\.'),'[',1)
dat <- subset(dat,!duplicated(paste(Normal,chromosome,start,REF,ALT)))
dat$refseq <- sapply(strsplit(dat$refseq,'\\.'),'[',1)
write.table(dat,file='lolliplot/RNF43_forLolliplot_0925.txt',row.names=F,quote=F,sep='\t')
by(dat$Normal,dat$Type,function(x){length(unique(x))})
dat$type <- ifelse(dat$class %in% c('Frame_Shift_Del','Frame_Shift_Ins'), 'Frameshift', ifelse(dat$class %in% c('Missense_Mutation'),'Missense', ifelse(dat$class %in% c('Nonsense_Mutation'),'Nonsense',ifelse(dat$class %in% c('Splice_Site'),'Splice site',NA))))
mutc <- as.matrix ( table(dat$type,dat$Note) )
mutc
fisher.test(matrix(c(13,2,9,4),nrow=2))
fisher.test(matrix(c(14,1,4,9),nrow=2))
> fisher.test(matrix(c(14,1,4,9),nrow=2))
        Fisher''s Exact Test for Count Data
data:  matrix(c(14, 1, 4, 9), nrow = 2)
p-value = 0.001068
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
    2.549145 1459.373323
sample estimates:
odds ratio
  26.71767
library(ggplot2)
a <- data.frame(type = rownames(mutc), sam = 'Early stage\nCRC', num = mutc[,1])
b <- data.frame(type = rownames(mutc), sam = 'mCRC', num = mutc[,2])
data <- rbind(a,b)
data$type <- factor(data$type, levels=rev(c('Frameshift','Nonsense','Splice site','Missense')), labels=rev(c('Frameshift','Nonsense','Splice site','Missense'))	)
p <- ggplot(data = data, aes(x = sam, y = num, fill = factor(type)))+
     geom_bar(stat="identity",position='fill',width=0.7) +
	 theme_classic() +
     scale_fill_manual(values = rev(c("#DB3C3C", "#FF7E0E", "#6633FF","#3986CC"))) +
	 xlab('') +
	 ylab(expression(paste("% of ", italic("RNF43"), " mutations"))) +
	 theme(legend.title=element_blank(),
	 axis.text.x = element_text(color='black'),
     axis.text.y = element_text(color='black'),
     axis.ticks = element_line(color='black'))
p
ggsave('RNF43/RNF43_mutType_plot_0925.pdf',p, width=3.2, height=2.2)
extracellular (EC) ; 1 87 ;
protease-associated (PA) ; 88 186 ;
transmembrane (TM) ; 187 272 ;
RING finger (RING) ; 273 313 ;
DVL2-binding (DVL) ; 314 596 ;
C-terminal ; 597 783 ;
a <- data.frame(sam='mCRC',
                domain = c('extracellular (EC)','protease-associated (PA)','transmembrane (TM)','RING finger (RING)','DVL2-binding (DVL)','C-terminal'),
                bpstart = c(1,88,187,273,314,597),
				bpend = c(87,186,272,313,596,783),
				num = c(6,3,5,0,0,1))
b <- data.frame(sam='Early stage\nCRC',
                domain = c('extracellular (EC)','protease-associated (PA)','transmembrane (TM)','RING finger (RING)','DVL2-binding (DVL)','C-terminal'),
                bpstart = c(1,88,187,273,314,597),
				bpend = c(87,186,272,313,596,783),
				num = c(4,0,0,2,4,3))
data <- rbind(a,b)
data$domain <- factor(data$domain, levels=rev(a$domain), labels=rev(a$domain))
p <- ggplot(data = data, aes(x = sam, y = num, fill = factor(domain)))+
     geom_bar(stat="identity",position='fill',width=0.7) +
	 theme_classic() +
     scale_fill_manual(values = rev(c("#9DC3E6", "#FFE288", "#8DD3C3", "#F09E64", "#A9D18E", "#8FAADC"))) +
	 xlab('') +
	 ylab(expression(paste("% of ", italic("RNF43"), " mutations"))) +
	 theme(legend.title=element_blank(),
	 axis.text.x = element_text(color='black'),
     axis.text.y = element_text(color='black'),
     axis.ticks = element_line(color='black'))
p
ggsave('RNF43/RNF43_mut_distribution_plot_0925.pdf',p, width=4.1, height=2.2)
