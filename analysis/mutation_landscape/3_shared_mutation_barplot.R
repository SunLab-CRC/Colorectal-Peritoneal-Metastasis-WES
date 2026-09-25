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
re <- data.frame()
for(i in unique(info$Normal)){
	sinf <- subset(info, Normal == i & Class != 'Primary')
	sda <- subset(dat, Normal.x == i)
	pri <- subset(sda, Class=='Primary')
	mut_id_pri <- unique(pri$mutation_id.x)
	mutn_pri <- length(mut_id_pri)
	abd <- subset(sda, Class=='Nodule')
	mut_id_abd <- unique(abd$mutation_id.x)
	mutn_abd <- length(mut_id_abd)
	shared <- unique(intersect(mut_id_pri, mut_id_abd))
	sre <- data.frame(Normal=i, PM = 'abdomen',
				primary.private = mutn_pri - length(shared),
				shared = length(shared),
				metastasis.private = mutn_abd - length(shared))
	re <- rbind(re, sre)
	rm(shared,sre,abd,mut_id_abd,mutn_abd)
	ln <- subset(sda, Class=='Lymph')
	if(nrow(ln)>0){
	mut_id_ln <- unique(ln$mutation_id.x)
	mutn_ln <- length(mut_id_ln)
    shared <- unique(intersect(mut_id_pri, mut_id_ln))
	sre <- data.frame(Normal=i, PM = 'lymph',
				primary.private = mutn_pri - length(shared),
				shared = length(shared),
				metastasis.private = mutn_ln - length(shared))
	re <- rbind(re, sre)
	rm(shared,sre,mut_id_ln,mutn_ln,ln)
	}
	rm(sinf,sda,pri,i,mut_id_pri, mutn_pri)
	}
write.csv(re,file='ALL_mutation_type.csv',row.names=F,quote=F)
re <- read.csv('ALL_mutation_type.csv')
dat <- merge(re,subset(info,!duplicated(Normal)),by='Normal')
dat$sam_name <- sapply(strsplit(dat$sam_name,' '),'[',1)
dat$met <- ifelse(dat$PM=='abdomen','A','LN')
dat$sam_name <- paste(dat$sam_name,dat$met)
dat$share.prop <- dat$shared/(dat$shared + dat$primary.private + dat$metastasis.private)
dat$met <- factor(dat$met,levels=c('LN','A'))
dat2 <- dat[order(dat$met,dat$share.prop,decreasing=T),]
dat1 <- rbind(data.frame(id=dat$sam_name, mutn=dat$primary.private, note='Primary tumor private'),
              data.frame(id=dat$sam_name, mutn=dat$shared, note='Shared'),
              data.frame(id=dat$sam_name, mutn=dat$metastasis.private, note='Metastasis private')
			  )
dat1$sam <- factor(dat1$id, levels=dat2$sam_name,labels=dat2$sam_name)
dat1$note <- factor(dat1$note, levels=c('Shared','Primary tumor private','Metastasis private'),labels=c('Shared','Primary tumor private','Metastasis private'))
p <- ggplot(data = dat1, aes(x = sam, y = mutn, fill = factor(note)))+
     geom_bar(stat="identity",position='fill',width=0.7) +
	 theme_classic() +
     scale_fill_manual(values = c("#C2A5CE", "#A6D49F","#762A82")) +
	 xlab('') + ylab('Percentage') +
	 theme(legend.position="none",
	 axis.text.x = element_text(color='black',angle = 60, hjust = 1, vjust = 1),
     axis.text.y = element_text(color='black'),
     axis.ticks = element_line(color='black'))
ggsave('NJMU_mCRC_ALL_mutation_PM_barplot.pdf',p, width=3.5, height=2.5)
a <- subset(dat2, PM!='lymph')
dat2$met.prop <- dat2$metastasis.private / (dat2$primary.private + dat2$shared + dat2$metastasis.private)
subset(dat2, sam_name %in% c('CASE_ID LN','CASE_ID A'))
subset(dat2, sam_name %in% c('CASE_ID LN','CASE_ID A'))
subset(dat2, sam_name %in% c('CASE_ID LN','CASE_ID A'))
col = c(rgb(red=50,green=123,blue=186,alpha=255,max=255),
		rgb(red=222,green=10,blue=21,alpha=255,max=255),
		rgb(red=241,green=125,blue=0,alpha=255,max=255),
		rgb(red=179,green=213,blue=230,alpha=255,max=255),
		rgb(red=69,green=47,blue=140,alpha=255,max=255),
		rgb(red=147,green=12,blue=19,alpha=255,max=255),
		rgb(red=65,green=174,blue=119,alpha=255,max=255))
colors = structure(col, names = c('Missense','Nonsense','Frameshift','In_frame_indel','Splice_site','Nonstop','Multiple_Hits'))
ldg_hp = Legend(labels = names(colors) , title = "Mutation" ,
	 legend_gp = gpar(fill = colors , bar_width = 1 , fontsize = 12) , ncol = 3 ,
	 gap = unit(0.6, "cm"),
	 row_gap = unit(2, "mm")
	)
col <- rev(c("#762A82", "#C2A5CE", "#A6D49F"))
col_mut <- structure(col, names = c('Primary tumor private','Shared','Metastasis private'))
ldg_mut = Legend(labels = names(col_mut) , title = "Mutation type" ,
	 legend_gp = gpar(fill = col_mut , bar_width = 1 , fontsize = 12) , ncol = 2 ,
	 gap = unit(0.6, "cm") ,
	 row_gap = unit(2, "mm")
	)
pdf('legend.pdf',height=6.9,width=4.8)
draw(ldg_hp,, x = unit(0.2, "npc"), y = unit(0.6, "npc"),just = c("left", "bottom"))
draw(ldg_mut,, x = unit(0.2, "npc"), y = unit(0.48, "npc"),just = c("left", "bottom"))
dev.off()
