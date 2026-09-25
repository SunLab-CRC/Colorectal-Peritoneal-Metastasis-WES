setwd('<PROJECT_ROOT>')
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(circlize)
rm(list=ls())
info <- read.csv('<PROJECT_ROOT>')
info$indsam <- paste(info$Normal, info$Class_a, sep=':')
info$Class <- sapply(strsplit(info$Class_a,'-'),'[',1)
maf <- fread('NJMU_mCRC_ALL.Gene_mutation.maf',data.table=F)
maf$indsam <- paste(maf$Normal, maf$Class_a, sep=':')
dat <- merge(maf, info,by='indsam')
dat <- subset(dat, variant_allele_frequency>0)
dat$Class <- sapply(strsplit(dat$Class_a.x,'-'),'[',1)
rm(maf)
dat <- subset(dat, Variant_Classification %in% c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site'))
dri <- read.csv('../Colorectal_Cancer_Mut_Driver_Genes.csv')
re <- data.frame()
for(i in unique(info$Normal)){
	sinf <- subset(info, Normal == i & Class != 'Primary')
	sda <- subset(dat, Normal.x == i)
	pri <- subset(sda, Class=='Primary')
	mut_pri <- data.frame(id = unique(pri$mutation_id.x))
	mut_pri$SYMBOL <- sapply(strsplit(mut_pri$id,':'),'[',1)
	mutn_pri <- length(unique(mut_pri$id))
	abd <- subset(sda, Class=='Nodule')
	mut_abd <- data.frame(id = unique(abd$mutation_id.x))
	mut_abd$SYMBOL <- sapply(strsplit(mut_abd$id,':'),'[',1)
	mutn_abd <- length(unique(mut_abd$id))
	shared <- subset(mut_abd, id %in% mut_pri$id)
	pp <- subset(mut_pri, !id %in% shared$id)
	ap <- subset(mut_abd, !id %in% shared$id)
	sre <- data.frame(Normal=i, PM = 'abdomen',
				num.primary.private = nrow(pp),
				num.driver.primary.private = nrow(subset(pp, SYMBOL %in% dri$SYMBOL)),
				num.shared = nrow(shared),
				num.driver.shared = nrow(subset(shared, SYMBOL %in% dri$SYMBOL)),
				num.metastasis.private = nrow(ap),
				num.driver.metastasis.private = nrow(subset(ap, SYMBOL %in% dri$SYMBOL)))
	re <- rbind(re, sre)
	rm(abd, mut_abd, mutn_abd, shared, pp, ap, sre)
	ln <- subset(sda, Class=='Lymph')
	if(nrow(ln)>0){
	mut_ln <- data.frame(id = unique(ln$mutation_id.x))
	mut_ln$SYMBOL <- sapply(strsplit(mut_ln$id,':'),'[',1)
	mutn_ln <- length(unique(mut_ln$id))
	shared <- subset(mut_ln, id %in% mut_pri$id)
	pp <- subset(mut_pri, !id %in% shared$id)
	lnp <- subset(mut_ln, !id %in% shared$id)
	sre <- data.frame(Normal=i, PM = 'lymph',
				num.primary.private = nrow(pp),
				num.driver.primary.private = nrow(subset(pp, SYMBOL %in% dri$SYMBOL)),
				num.shared = nrow(shared),
				num.driver.shared = nrow(subset(shared, SYMBOL %in% dri$SYMBOL)),
				num.metastasis.private = nrow(lnp),
				num.driver.metastasis.private = nrow(subset(lnp, SYMBOL %in% dri$SYMBOL)))
	re <- rbind(re, sre)
	rm(ln, mut_ln, mutn_ln, shared, pp, lnp, sre)
	}
	rm(sinf,sda,pri,i,mut_pri, mutn_pri)
	}
write.csv(re,file='ALL_mutation_Driver_enrich_data.csv',row.names=F,quote=F)
re <- read.csv('ALL_mutation_Driver_enrich_data.csv')
re$pri.pri.prop <- re$num.driver.primary.private / re$num.primary.private
re$shared.prop <- re$num.driver.shared / re$num.shared
re$met.prop <- re$num.driver.metastasis.private / re$num.metastasis.private
gencode <- read.table('gencode.v19.annotation.bed')
re$pri.enrich <- re$pri.pri.prop / (139/20345)
re$shared.enrich <- re$shared.prop / (139/20345)
re$met.enrich <- re$met.prop / (139/20345)
wilcox.test(re$pri.enrich, re$shared.enrich, paired = TRUE)
wilcox.test(re$met.enrich, re$shared.enrich, paired = TRUE)
wilcox.test(re$met.enrich, re$pri.enrich, paired = TRUE)
dat <- subset(re, PM == 'abdomen')
wilcox.test(dat$pri.enrich, dat$shared.enrich, paired = TRUE)
wilcox.test(dat$met.enrich, dat$shared.enrich, paired = TRUE)
wilcox.test(dat$met.enrich, dat$pri.enrich, paired = TRUE)
library(ggplot2)
library(ggprism)
dat1 <- rbind(data.frame(id=dat$Normal, mutn=dat$pri.enrich, note='Primary\ntumor private'),
              data.frame(id=dat$Normal, mutn=dat$shared.enrich, note='Shared'),
              data.frame(id=dat$Normal, mutn=dat$met.enrich, note='Metastasis\nprivate')
			  )
dat1$note <- factor(dat1$note, levels=c('Shared','Primary\ntumor private','Metastasis\nprivate'), labels=c('Shared','Primary\ntumor private','Metastasis\nprivate'))
attach(dat1)
p <- ggplot(data = dat1, aes(x = note, y = mutn))+
     geom_violin(width=0.9,aes(fill = note)) +
	 geom_boxplot(width=0.1, color="black", alpha=0) +
	 theme_classic() +
     scale_fill_manual(values = c("#C2A5CE", "#A6D49F","#762A82")) +
	 xlab('') + ylab('Fold enrichment in\nCRC driver genes') +
	 ylim(0,25) +
	 theme(legend.position="none",
	 axis.text.x = element_text(color='black',angle = 30, hjust = 1, vjust = 1),
     axis.text.y = element_text(color='black'),
     axis.ticks = element_line(color='black'))
p_vals <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~label,
  "Shared",   "Primary\ntumor private", 0.008, 24, "p.exprs.ital",
  "Primary\ntumor private",   "Metastasis\nprivate", 1, 17, "p.exprs.ital"
)
p1 = p + add_pvalue(p_vals, label = "p-value = {p.adj}", tip.length = 0.02, label.size = 3)
p1
ggsave('NJMU_mCRC_ALL_mutation_Driver.enrich_violinPlot.pdf',p1, width=2.8, height=2.5)
