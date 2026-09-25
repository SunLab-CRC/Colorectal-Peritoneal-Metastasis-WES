options(stringsAsFactors=F)
library(data.table)
library(ggplot2)
setwd('<PROJECT_ROOT>')
rm(list=ls())
pri <- read.csv('MSK_Abdomen_Metastasis_SCNA.Amplification.Driver_Gscore.csv')
met <- read.csv('MSK_TCGA_Early_Primary_SCNA.Amplification.Driver_Gscore.csv')
dat <- merge(met[,c('ID','symbol','chr.amp','start.amp','end.amp','log10Q.met','G.score.met')], pri[,c('ID','log10Q.pri','G.score.pri')], by='ID')
cgen <- read.csv('../../2_driverGene_waterfall_plot/SCNA_driver_CRC_Amp_candidate_Gene.csv')
cgen$ID <- sapply(strsplit(cgen$ID,'\\.'),'[',1)
cgen$Symbol <- toupper(cgen$Symbol)
dat$Symbol <- toupper(dat$symbol)
dat <- subset(dat, symbol %in% cgen$Symbol | ID %in% cgen$ID)
dat <- subset(dat,!duplicated(ID))
dat$Gscore_diff <- dat$G.score.met - dat$G.score.pri
sig <- subset(dat, log10Q.met > -log10(0.05) )
dat$note <- ifelse(dat$log10Q.met > -log10(0.05), 'sig', 'non.sig')
dat$log10Q.met <- ifelse(dat$log10Q.met==Inf, 25, dat$log10Q.met)
library(ggrepel)
sig$symbol <- paste("italic('",sig$symbol,"')",sep='')
p <- ggplot(data = dat, aes(x = Gscore_diff, y = log10Q.met, color=note))+
     geom_point(size=0.5) +
	 theme_classic()  +
	 scale_color_manual(values = rev(c("#CC0000", "#000000"))) +
	 xlab('GISTIC score difference') +
	 ylab(expression(paste(-log[10]," (FDR)"))) +
	 geom_hline(yintercept=-log10(0.05),color='grey50',linetype="dashed") +
	 theme(legend.position="none",
	 axis.text.x = element_text(color='black'),
     axis.text.y = element_text(color='black'),
     axis.ticks = element_line(color='black')) +
	 xlim(-0.11,0.15) + ylim(0,28) +
	 geom_text_repel(data = sig,
	 parse = TRUE,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                  aes(label = symbol),
				  color = '#CC0000',
				  point.padding = 0.3,
                  size = 3) +
				  labs(title = "Amplification",size=2)
ggsave('MSK_TCGA_Validation_SCNA.AMP.Driver.pdf',p,width=4.5, height=3.5)
options(stringsAsFactors=F)
library(data.table)
library(ggplot2)
setwd('<PROJECT_ROOT>')
rm(list=ls())
pri <- read.csv('MSK_TCGA_Early_Primary_SCNA.Deletion.Driver_Gscore.csv')
met <- read.csv('MSK_Abdomen_Metastasis_SCNA.Deletion.Driver_Gscore.csv')
dat <- merge(met[,c('ID','symbol','chr.del','start.del','end.del','log10Q.met','G.score.met')], pri[,c('ID','log10Q.pri','G.score.pri')], by='ID')
cgen <- read.csv('../../2_driverGene_waterfall_plot/SCNA_driver_CRC_Del_candidate_Gene.csv')
cgen$ID <- sapply(strsplit(cgen$ID,'\\.'),'[',1)
cgen$Symbol <- toupper(cgen$Symbol)
dat$Symbol <- toupper(dat$symbol)
dat <- subset(dat, symbol %in% cgen$Symbol | ID %in% cgen$ID)
dat <- subset(dat,!duplicated(ID))
dat$Gscore_diff <- dat$G.score.met - dat$G.score.pri
sig <- subset(dat, log10Q.met > -log10(0.05) )
dat$note <- ifelse(dat$log10Q.met > -log10(0.05), 'sig', 'non.sig')
dat$log10Q.met <- ifelse(dat$log10Q.met==Inf, 20, dat$log10Q.met)
library(ggrepel)
sig$symbol <- paste("italic('",sig$symbol,"')",sep='')
p <- ggplot(data = dat, aes(x = Gscore_diff, y = log10Q.met, color=note))+
     geom_point(size=0.5) +
	 theme_classic()  +
	 scale_color_manual(values = rev(c("#000080", "#000000"))) +
	 xlab('GISTIC score difference') +
	 ylab(expression(paste(-log[10]," (FDR)"))) +
	 geom_hline(yintercept=-log10(0.05),color='grey50',linetype="dashed") +
	 theme(legend.position="none",
	 axis.text.x = element_text(color='black'),
     axis.text.y = element_text(color='black'),
     axis.ticks = element_line(color='black')) +
	 xlim(-0.2,0.16) + ylim(0,20.5) +
	 geom_text_repel(data = sig,
	 parse = TRUE,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                  aes(label = symbol),
				  color = '#000080',
				  point.padding = 0.3,
                  size = 3) +
				  labs(title = "Losses",size=2)
ggsave('MSK_TCGA_Validation_SCNA.DEL.Driver.pdf',p,width=4.5, height=3.5)
p <- ggplot(data = dat, aes(x = Gscore_diff, y = log10Q.met, color=note))+
     geom_point(size=0.5) +
	 theme_classic()  +
	 scale_color_manual(values = rev(c("#0066CC", "#000000"))) +
	 xlab('GISTIC score difference') + ylab('-log10(FDR)') +
	 geom_hline(yintercept=-log10(0.05),color='grey50',linetype="dashed") +
	 theme(legend.position="none",
	 axis.text.x = element_text(color='black'),
     axis.text.y = element_text(color='black'),
     axis.ticks = element_line(color='black')) +
	 xlim(-0.25,0.25) +
	  geom_text_repel(data = subset(dat, !symbol %in% c('DCC','FBXW7','FGFR3','MAP2K4','SMAD4','CCNC') & log10Q.met > -log10(0.05) ),
	 min.segment.length = Inf,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = symbol),
				  color = 'grey50',
				  point.padding = 0.3,
                  size = 2) +
	 geom_text_repel(data = subset(dat, symbol %in% c('DCC','FBXW7','FGFR3','MAP2K4','SMAD4','CCNC') & log10Q.met > -log10(0.05) ),
	 min.segment.length = Inf,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                  aes(label = symbol),
				  color = '#0066CC',
				  point.padding = 0.3,
                  size = 2) +
				  labs(title = "Deletion")
ggsave('MSK_TCGA_Validation_SCNA.DEL.Driver.pdf',p,width=3, height=2.7)
