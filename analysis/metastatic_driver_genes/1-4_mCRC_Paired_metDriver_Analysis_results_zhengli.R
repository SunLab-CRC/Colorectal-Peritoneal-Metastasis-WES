setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
maf <- read.csv("NJMU_mCRC_mutation_clonality_primary_metastases_vaf005_read4_0.5_0.1_corrected.csv")
maf$gene <- sapply(strsplit(maf$mutation_id,':'),'[',1)
maf$Variant_Classification <- sapply(strsplit(maf$mutation_id,':'),'[',2)
maf <- subset(maf, Variant_Classification %in% c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site'))
id <- read.csv('<PROJECT_ROOT>')
id$sam <- sapply(strsplit(id$sam_name,' '),'[',1)
id <- subset(id, !duplicated(sam))
mut <- merge(maf,id[,c('sam','Normal')],by='Normal')
gen_name = 'SMARCB1'
subset(mut,gene==gen_name)
data <- read.csv('NJMU_mCRC_Metastasis_Driver_Gene_Paired_Analysis_Results_vaf005_0.5_0.1.csv')
data$multinomial.chisq.FDR <- p.adjust(data$multinomial.chisq, method='BH')
subset(data, multinomial.chisq.FDR<0.05)
data <- subset(data, multinomial.chisq<0.05)
library(ggplot2)
data <- data[order((data$metastasis.prop), (1-data$primary.prop), data$no.mut, decreasing=T),]
data$gene = paste(data$gene,' (',round(data$multinomial.chisq,3),', ',data$no.mut,')',sep='')
dat1 <- rbind(data.frame(gene=data$gene, prop=data$metastasis.prop, note='Metastasis'),
              data.frame(gene=data$gene, prop=data$maintained.prop
, note='maintained'),
              data.frame(gene=data$gene, prop=data$primary.prop
, note='primary'),
              data.frame(gene=rep('Background',3),prop=c(0.3191863,0.3103968,0.3704169),note=c('maintained','Metastasis','primary'))
			  )
dat1$gene <- factor(dat1$gene, levels=rev(c('Background',data$gene)))
dat1$note <- factor(dat1$note, levels=rev(c('Metastasis','maintained','primary')), labels=rev(c('Metastasis favoured','Maintained','Primary favoured')))
ggplot(data = dat1, aes(x = gene, y = prop, fill = factor(note)))+
     geom_bar(stat="identity",position='fill',width=0.8) +
	 theme_classic() +
	 coord_flip() +
     scale_fill_manual(values = rev(c("#762A82", "#C2A5CE", "#A6D49F"))) +
	 xlab('') + ylab('') +
	 theme(legend.title=element_blank(),
	 axis.text.x = element_text(color='black'),
     axis.text.y = element_text(color='black'),
     axis.ticks = element_line(color='black'))
p <- ggplot(data = dat1, aes(x = gene, y = prop, fill = factor(note)))+
     geom_bar(stat="identity",position='fill',width=0.8) +
	 theme_classic() +
	 coord_flip() +
     scale_fill_manual(values = rev(c("#762A82", "#C2A5CE", "#A6D49F"))) +
	 xlab('') + ylab('Proportion of\nmutations') +
	 theme(legend.position="none",
	 axis.text.x = element_text(color='black'),
     axis.ticks = element_line(color='black'),
	 axis.line.y=element_blank(),
     axis.text.y=element_blank(),
     axis.ticks.y=element_blank(),
     axis.title.y=element_blank())
dat2 <- data[,c('gene','multinomial.chisq','no.mut')]
dat2$gene <- sapply(strsplit(dat2$gene,' '),'[',1)
dat2 <- rbind(data.frame(gene=NA,multinomial.chisq=NA, no.mut=NA),dat2)
dat2$sequ <- rev(1:nrow(dat2))
dat2$x <- 1
dat2$note <- paste('(',dat2$no.mut,')',sep='')
dat2$note[1] <- ''
dat2$gene[1] <- 'Background'
dat2$gene <- factor(dat2$gene, levels=rev(dat2$gene), labels=rev(dat2$gene))
mutn <- ggplot(data = dat2, aes(x = x, y = gene))+
     geom_text(label=dat2$note,hjust = 0) +
	 theme_classic() +
	 xlim(0.8,1.2) +
	 xlab('') + ylab('') +
	 theme(legend.position="none",
	 axis.line=element_blank(),
     axis.text=element_blank(),
     axis.ticks=element_blank(),
	 plot.margin = unit(c(0,0,0,0),"inches"))
pval <- ggplot(data = dat2, aes(x = -log10(multinomial.chisq), y = gene))+
     geom_segment(x=0,y=dat2$gene,xend=-log10(dat2$multinomial.chisq),yend=dat2$gene,color="#E0E0E0")  +
	 geom_point(size=2,color="#57B2AB") +
	 theme_classic() +
	 xlab(expression(paste(-log[10]," (", italic("P"), ")"))) +
	 ylab('') +
	 xlim(0,15) +
	 theme(legend.position="none",
	 axis.text = element_text(color='black'),
     axis.ticks = element_line(color='black'),
	 axis.ticks.y = element_blank(),
	 axis.text.y = element_text(face = "italic"))
library(patchwork)
complot <- pval + p + mutn + plot_layout(ncol = 3,width=c(2,2.2,1.2))
complot
ggsave('NJMU_mCRC_Paired_test_plot.pdf',complot, width=4.5, height=4)
