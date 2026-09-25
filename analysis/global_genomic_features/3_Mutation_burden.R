R
options(stringsAsFactors=F)
setwd('<PROJECT_ROOT>')
library(data.table)
library(maftools)
rm(list=ls())
work_dir <- "<PROJECT_ROOT>"
info <- fread(paste(work_dir,'config/tumor_normal.list',sep=''),data.table=F)
info <- subset(info, Normal != 'CASE_ID')
tre <- read.csv('<PROJECT_ROOT>')
tre$Class <- sapply(strsplit(tre$Class_a,'-'),'[',1)
tmb_qc <- data.frame()
for(n in unique(info$Normal)){
  sinf <- subset(info,Normal==n)
  stre <- subset(tre,Normal==n)
  sinf <- subset(sinf, Class_a %in% stre$Class_a)
  rm(stre)
  for(i in 1:nrow(sinf)){
    maf <- read.csv(paste(work_dir,'maf/',sinf$Tumor[i],"_",sinf$Normal[i], "_GGA_Filter_funcotated.maf",sep=''),comment.char = "#",sep='\t')
	maf$mutation_id <- paste(maf$Hugo_Symbol,maf$Variant_Classification,maf$Chromosome,maf$Start_Position,maf$Tumor_Seq_Allele1,maf$Tumor_Seq_Allele2,sep=':')
	cln <- read.delim(paste(work_dir,'Pyclone/result/Input_all_VAF005_read4/',n, "/tables/loci.tsv",sep=''),comment.char = "#",sep='\t')
    clu_use <- read.delim(paste('<PROJECT_ROOT>',n,'_cluster_use.txt',sep=''),h=F)
    cln <- merge(cln, clu_use,by='cluster_id')
	cln <- subset(cln, sample_id == sinf$Class_a[i] & variant_allele_frequency>0)
    smaf <- subset(maf, mutation_id %in% cln$mutation_id)
	write.table(smaf, file='tmp.maf',row.names=F,quote=F,sep='\t')
	smaf <- read.maf(maf = "tmp.maf")
	stmb <- as.data.frame(tmb(maf = smaf, logScale = F))
	system('rm tmp.maf')
	stmb$Normal <- sinf$Normal[i]
	stmb$Tumor <- sinf$Tumor[i]
	stmb$Class_a <- sinf$Class_a[i]
	tmb_qc <- rbind(tmb_qc, stmb)
	rm(i, maf, cln, clu_use, smaf, stmb)
	}
	rm(n, sinf)
  }
write.csv(tmb_qc, file="TMB_each_Sample_QC.csv",row.names=F,quote=F)
system('sz TMB_each_Sample_QC.csv')
options(stringsAsFactors=F)
setwd('<PROJECT_ROOT>')
library(data.table)
library(ggplot2)
library(RColorBrewer)
rm(list=ls())
sam <- read.csv('TMB_each_Sample_QC.csv')
id <- read.csv('../../0_patient_sample_info/Clonevol_sample_used.csv')
sam <- subset(sam, paste(Normal,Class) %in% paste(id$Normal,id$Class_a))
sam$Class_a <- sam$Class
sam$Class <- sapply(strsplit(sam$Class,'-'),'[',1)
dat <- data.frame()
for(i in unique(sam$Normal)){
	sams <- subset(sam, Normal == i)
	p <- subset(sams, Class == 'Primary')
	primary_mean <- mean(p$total_perMB)
	rm(p)
	n <- subset(sams, Class == 'Nodule')
	Nodule_mean <- mean(n$total_perMB)
	rm(n)
	l <- subset(sams, Class == 'Lymph')
	if(nrow(l)==0){
	  Lymph_mean <- NA
	}else{
	  Lymph_mean <- mean(l$total_perMB)
	}
	rm(l)
	m <- subset(sams, Class %in% c('Lymph','Nodule'))
	Met_mean <- mean(m$total_perMB)
	rm(m)
	dat <- rbind(dat,
	           data.frame(id = i,
			              primary_mean = primary_mean,
						  Nodule_mean  = Nodule_mean,
						  Lymph_mean   = Lymph_mean,
						  Met_mean     = Met_mean
						  ))
}
wilcox.test(dat$primary_mean, dat$Nodule_mean, paired = TRUE)
cor.test(dat$primary_mean, dat$Nodule_mean)
library(ggplot2)
library(ggprism)
setheme <- theme(axis.text.y = element_text(size=16, color = "black"),
		  axis.text.x = element_text(size=16, color = "black"),
		  axis.title.y = element_text(size=16, color = "black", face = "bold"),
		  axis.title.x = element_text(size=16, color = "black", face = "bold"),
		  plot.title = element_text(size=16, color = "black",hjust = 0.5,margin = margin(t = 1, r = 0, b = 0.5, l = 0,"cm"), , face = "bold")) +
	theme(axis.line = element_line(colour = 'black', linewidth = 0.5)) +
	theme(plot.margin = margin(t=0.2, r=2, b=0.2, l=1, "cm"))
p <- subset(sda,Class=='P')
n <- subset(sda,Class=='N')
pn <- merge(p,n,by='id')
p0 <- ggplot(aes(x = tmb.x, y = tmb.y),data=pn) +
      geom_point(fill='#3F76B4',color='#212121',alpha=0.6, size=4, shape=21) +
	  theme_classic() +
	  ylab("TMB in abdominal\nmetastases") +
	  xlab("TMB in primary CRC") +
	  ggtitle("   ") +
	  stat_smooth(method = "lm",formula = y ~ x,colour = "lightgrey",fill='lightgrey',alpha=0.6) +
	  setheme
p1 <- p0 +
      geom_text(x = 1.5, y = 3.2,label=expression(paste(r^2,"=0.82",sep='')), size=5) +
	  geom_text(x = 1.5, y = 2.9,label=expression(paste(italic(P),"=0.002",sep='')), size=5)
p1
ggsave('TMB_Primary_Nodule_scatterPlot.pdf',p1, width=5 ,height= 3.8)
sda <- rbind(data.frame(id=dat$id, rate=dat$primary_mean, Class='P'),
  data.frame(id=dat$id, rate=dat$Nodule_mean, Class='N'),
  data.frame(id=dat$id, rate=dat$Lymph_mean, Class='LN')
  )
sda <- subset(sda, !is.na(rate))
setheme <- theme(axis.text.y = element_text(size=16, color = "black"),
		  axis.text.x = element_text(size=16, color = "black"),
		  axis.title.y = element_text(size=16, color = "black", face = "bold"),
		  axis.title.x = element_text(size=16, color = "black", face = "bold"),
		  plot.title = element_text(size=16, color = "black",hjust = 0.5,margin = margin(t = 1, r = 0, b = 0.5, l = 0,"cm"), , face = "bold")) +
	theme(axis.line = element_line(colour = 'black', size = 0.5)) +
	theme(plot.margin = margin(t=0.2, r=2, b=0.2, l=1, "cm"))
sda$Type <- ifelse(sda$Class=='P','P\n(n=11)',ifelse(sda$Class=='LN','LN\n(n=3)',ifelse(sda$Class=='N','A\n(n=11)',NA)))
sda$Type <- factor(sda$Type,levels=c('P\n(n=11)','A\n(n=11)','LN\n(n=3)'),labels=c('P\n(n=11)','A\n(n=11)','LN\n(n=3)'))
p0 <- ggplot(aes(x=Type,y=rate),data=sda) +
      geom_boxplot(fill='#3F76B4',color='#212121',width = 0.6,alpha=0.2,outlier.shape = NA) +
	  geom_jitter(width=0.1,color='#212121',size=3.5,alpha=0.3) +
	  theme_classic() +
	  xlab('') +
	  ylab('TMB') +
	  ylim(0.8,4.5) +
	  setheme
wilcox.test(dat$primary_mean, dat$Nodule_mean, paired = TRUE)
wilcox.test(dat$primary_mean, dat$Lymph_mean, paired = TRUE)
p_vals <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~label,
  "P\n(n=11)",   "A\n(n=11)",     0.21, 3.8, "p.exprs.ital",
  "P\n(n=11)",   "LN\n(n=3)",     1.00, 4.3, "p.exprs.ital"
)
p1 = p0 + add_pvalue(p_vals, label = "P = {p.adj}", tip.length = 0.02, label.size = 5)
p1
ggsave('TMB_Primary.VS.Met_paired.pdf',p1, width=5 ,height= 3.6)
wilcox.test(total_perMB ~ Class, data = subset(sam, Class %in% c('Primary','Nodule')))
wilcox.test(total_perMB ~ Class, data = subset(sam, Class %in% c('Primary','Lymph')))
setheme <- theme(axis.text.y = element_text(size=16, color = "black"),
		  axis.text.x = element_text(size=16, color = "black"),
		  axis.title.y = element_text(size=16, color = "black", face = "bold"),
		  axis.title.x = element_text(size=16, color = "black", face = "bold"),
		  plot.title = element_text(size=16, color = "black",hjust = 0.5,margin = margin(t = 1, r = 0, b = 0.5, l = 0,"cm"), , face = "bold")) +
	theme(axis.line = element_line(colour = 'black', size = 0.5)) +
	theme(plot.margin = margin(t=0.2, r=2, b=0.2, l=1, "cm"))
sam$Type <- ifelse(sam$Class=='Primary','Primary CRC\n(n=16)',ifelse(sam$Class=='Lymph','Lymphatic metastases\n(n=4)',ifelse(sam$Class=='Nodule','Abdominal metastases\n(n=29)',NA)))
sam$Type <- factor(sam$Type,levels=c('Primary CRC\n(n=16)','Abdominal metastases\n(n=29)','Lymphatic metastases\n(n=4)'),labels=c('Primary\nCRC\n(n=16)','Abdominal\nmetastases\n(n=29)','Lymphatic\nmetastases\n(n=4)'))
p0 <- ggplot(aes(x=Type,y=total_perMB),data=sam) +
      geom_boxplot(fill='#3F76B4',color='#212121',width = 0.6,alpha=0.2,outlier.shape = NA) +
	  geom_jitter(width=0.1,color='#212121',size=3.5,alpha=0.3) +
	  theme_classic() +
	  xlab('') +
	  ylim(0,4.5) +
	  ylab('TMB') +
	  setheme
p_vals <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~label,
  "Primary\nCRC\n(n=16)",   "Abdominal\nmetastases\n(n=29)", 0.04, 3.9, "p.exprs.ital",
  "Primary\nCRC\n(n=16)",   "Lymphatic\nmetastases\n(n=4)", 0.32, 4.3, "p.exprs.ital"
)
p1 = p0 + add_pvalue(p_vals, label = "P = {p.adj}", tip.length = 0.02, label.size = 5)
p1
ggsave('TMB_Primary.VS.Nodule_region.pdf',p1, width=6 ,height= 3.8)
by(sam$total_perMB, sam$Class, median)
