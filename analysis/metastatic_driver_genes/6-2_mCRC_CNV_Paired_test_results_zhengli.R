setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
library(data.table)
library(ggplot2)
rm(list=ls())
dat <- read.csv('NJMU_Amplification.Metastasis.Driver.Gene_Paired_Analysis_Results_cytoband.csv')
dat$no.mut <- dat$maintained + dat$metastasis.favored + dat$primary.favored
dat <- subset(dat, no.mut>=5)
dat$metastasis.prop <- dat$metastasis.favored / dat$no.mut
dat$maintained.prop <- dat$maintained / dat$no.mut
dat$primary.prop <- dat$primary.favored / dat$no.mut
dat$logFC.metastasis <- log2(dat$metastasis.prop / 0.2741716 )
dat$logFC.maintained <- log2(dat$maintained.prop / 0.3118664)
dat$logFC.primary <- log2(dat$primary.prop / 0.4139620)
dat$multinomial.chisq.FDR <- p.adjust(dat$multinomial.chisq, method='BH')
amp <- subset(dat, multinomial.chisq < 0.05 & logFC.primary<0)
cgen <- read.csv('../../2_driverGene_waterfall_plot/SCNA_driver_CRC_Amp_candidate_Gene.csv')
cgen$ID <- sapply(strsplit(cgen$ID,'\\.'),'[',1)
cgen$Symbol <- toupper(cgen$Symbol)
amp$ID <- sapply(strsplit(amp$ID,'\\.'),'[',1)
amp$Symbol <- toupper(sapply(strsplit(amp$gene,' '),'[',2))
amp <- subset(amp, Symbol %in% cgen$Symbol | ID %in% cgen$ID)
amp <- subset(amp,!duplicated(ID))
write.csv(amp,file='NJMU_Metastasis.Driver.Amplification_Paired_Analysis_sigResults.csv',row.names=F,quote=F)
rm(list=ls())
library(ggplot2)
data <- read.csv('NJMU_Metastasis.Driver.Amplification_Paired_Analysis_sigResults.csv')
data <- data[order((data$metastasis.prop), (1-data$primary.prop), data$no.mut, decreasing=T),]
data$gene <- data$Symbol
dat1 <- rbind(data.frame(gene=data$gene, prop=data$metastasis.prop, note='Metastasis'),
              data.frame(gene=data$gene, prop=data$maintained.prop
, note='maintained'),
              data.frame(gene=data$gene, prop=data$primary.prop
, note='primary'),
              data.frame(gene=rep('Background',3),prop=c(0.3118664, 0.2741716,0.4139620),note=c('maintained','Metastasis','primary'))
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
	 xlab('') + ylab('') +
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
mutn <- ggplot(data = dat2, aes(x = x, y = gene))+
     geom_text(label=dat2$note,hjust = 0) +
	 theme_classic() +
	 xlim(0.8,1.2) +
	 xlab('') + ylab('') +
	 theme(legend.position="none",
	 axis.line=element_blank(),
     axis.text=element_blank(),
     axis.ticks=element_blank())
dat2$gene[1] <- 'Background'
dat2$gene <- factor(dat2$gene, levels=rev(dat2$gene), labels=rev(dat2$gene))
pval <- ggplot(data = dat2, aes(x = -log10(multinomial.chisq), y = gene))+
     geom_segment(x=0,y=dat2$gene,xend=-log10(dat2$multinomial.chisq),yend=dat2$gene,color="#E0E0E0")  +
	 geom_point(size=2,color="#57B2AB") +
	 theme_classic() +
	 xlab(expression(paste(-log[10]," (", italic("P"), ")"))) +
	 ylab('') +
	 theme(legend.position="none",
	 axis.text = element_text(color='black'),
     axis.ticks = element_line(color='black'),
	 axis.ticks.y = element_blank(),
	 axis.text.y = element_text(face = "italic"))
library(patchwork)
complot <- pval + p + mutn + plot_layout(ncol = 3)
ggsave('NJMU_mCRC_sCNA.Amplification_Paired_test_plot.pdf',complot, width=6, height=2.5)
rm(list=ls())
dat <- read.csv('NJMU_Deletion.Metastasis.Driver.Gene_Paired_Analysis_Results_cytoband.csv')
dat$no.mut <- dat$maintained + dat$metastasis.favored + dat$primary.favored
dat <- subset(dat, no.mut>=5)
dat$metastasis.prop <- dat$metastasis.favored / dat$no.mut
dat$maintained.prop <- dat$maintained / dat$no.mut
dat$primary.prop <- dat$primary.favored / dat$no.mut
dat$logFC.metastasis <- log2(dat$metastasis.prop / 0.2849778 )
dat$logFC.maintained <- log2(dat$maintained.prop / 0.5556084)
dat$logFC.primary <- log2(dat$primary.prop / 0.1594138)
dat$multinomial.chisq.FDR <- p.adjust(dat$multinomial.chisq, method='BH')
del <- subset(dat, multinomial.chisq < 0.05 & logFC.primary<0)
cgen <- read.csv('../../2_driverGene_waterfall_plot/SCNA_driver_CRC_Del_candidate_Gene.csv')
cgen$ID <- sapply(strsplit(cgen$ID,'\\.'),'[',1)
cgen$Symbol <- toupper(cgen$Symbol)
del$ID <- sapply(strsplit(del$ID,'\\.'),'[',1)
del$Symbol <- toupper(sapply(strsplit(del$gene,' '),'[',2))
del <- subset(del, Symbol %in% cgen$Symbol | ID %in% cgen$ID)
del <- subset(del,!duplicated(ID))
write.csv(del,file='NJMU_Metastasis.Driver.Deletion_Paired_Analysis_sigResults.csv',row.names=F,quote=F)
rm(list=ls())
library(ggplot2)
data <- read.csv('NJMU_Metastasis.Driver.Deletion_Paired_Analysis_sigResults.csv')
data <- data[order((data$metastasis.prop), (1-data$primary.prop), data$no.mut, decreasing=T),]
data$gene <- data$Symbol
dat1 <- rbind(data.frame(gene=data$gene, prop=data$metastasis.prop, note='Metastasis'),
              data.frame(gene=data$gene, prop=data$maintained.prop
, note='maintained'),
              data.frame(gene=data$gene, prop=data$primary.prop
, note='primary'),
              data.frame(gene=rep('Background',3),prop=c(0.5556084,0.2849778,0.1594138),note=c('maintained','Metastasis','primary'))
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
     axis.ticks = element_line(color='black'),
	 axis.text.y = element_text(color='black',face = "italic"))
p <- ggplot(data = dat1, aes(x = gene, y = prop, fill = factor(note)))+
     geom_bar(stat="identity",position='fill',width=0.8) +
	 theme_classic() +
	 coord_flip() +
     scale_fill_manual(values = rev(c("#762A82", "#C2A5CE", "#A6D49F"))) +
	 xlab('') + ylab('') +
	 theme(legend.position="none",
	 axis.text.x = element_text(color='black'),
     axis.ticks = element_line(color='black'),
	 axis.line.y=element_blank(),
     axis.ticks.y=element_blank(),
     axis.title.y=element_blank(),
	 axis.text.y = element_blank())
dat2 <- data[,c('gene','multinomial.chisq','no.mut')]
dat2$gene <- sapply(strsplit(dat2$gene,' '),'[',1)
dat2 <- rbind(data.frame(gene=NA,multinomial.chisq=NA, no.mut=NA),dat2)
dat2$sequ <- rev(1:nrow(dat2))
dat2$x <- 1
dat2$note <- paste('(',dat2$no.mut,')',sep='')
dat2$note[1] <- ''
mutn <- ggplot(data = dat2, aes(x = x, y = gene))+
     geom_text(label=dat2$note,hjust = 0) +
	 theme_classic() +
	 xlim(0.8,1.2) +
	 xlab('') + ylab('') +
	 theme(legend.position="none",
	 axis.line=element_blank(),
     axis.text=element_blank(),
     axis.ticks=element_blank())
dat2$gene[1] <- 'Background'
dat2$gene <- factor(dat2$gene, levels=rev(dat2$gene), labels=rev(dat2$gene))
pval <- ggplot(data = dat2, aes(x = -log10(multinomial.chisq), y = gene))+
     geom_segment(x=0,y=dat2$gene,xend=-log10(dat2$multinomial.chisq),yend=dat2$gene,color="#E0E0E0")  +
	 geom_point(size=2,color="#57B2AB") +
	 theme_classic() +
	 xlab(expression(paste(-log[10]," (", italic("P"), ")"))) +
	 ylab('') +
	 theme(legend.position="none",
	 axis.text = element_text(color='black'),
     axis.ticks = element_line(color='black'),
	 axis.ticks.y = element_blank(),
	 axis.text.y = element_text(face = "italic"))
library(patchwork)
complot <- pval + p + mutn + plot_layout(ncol = 3)
ggsave('NJMU_mCRC_sCNA.Deletion_Paired_test_plot.pdf',complot, width=6, height=2.5)
