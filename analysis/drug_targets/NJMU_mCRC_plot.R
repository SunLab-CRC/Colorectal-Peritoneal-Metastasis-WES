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
snv <- read.csv('NJMU_mCRC_OncoKB_mutation_forPlot.csv')
snv <- snv[,1:9]
snv$indsam <- gsub('_','-',paste(snv$Normal, snv$Class_a, sep=':'))
snv <- merge(snv, info[,c('indsam','sam_name')], by='indsam')
cnv <- read.csv('NJMU_mCRC_OncoKB_sCNA_forPlot.csv')
cnv <- subset(cnv, nchar(id)>0)
cnv <- cnv[,1:6]
sam <- read.csv('<PROJECT_ROOT>')
sam$id <- paste(sam$Tumor,sam$Normal,sep='_')
cnv <- merge(cnv, sam[,c('id','Normal','Class_a')], by='id')
cnv$indsam <- paste(cnv$Normal, cnv$Class_a, sep=':')
subset(cnv, !indsam %in% info$indsam)
cnv <- merge(cnv, info[,c('indsam','sam_name')], by='indsam')
snv <- snv[,c('Hugo_Symbol','Variant_Classification','level','sam_name')]
cnv <- cnv[,c('gene','Alterations','level','sam_name')]
intersect(snv$Hugo_Symbol, cnv$Hugo_Symbol)
dat <- rbind(snv,cnv)
dat$sam <- sapply(strsplit(dat$sam_name,' '),'[',1)
a <- subset(dat, level != 'R1')
a <- a[order(a$Hugo_Symbol, a$sam_name, a$Variant_Classification, a$level),]
a <- subset(a, !duplicated(paste(a$Hugo_Symbol, a$Variant_Classification, a$sam_name)))
r <- subset(dat, level == 'R1')
data <- rbind(a,r)
a <- a[order(a$sam, a$level),]
sam_dt <- subset(a, !duplicated(sam))
re <- data.frame(drug = c('Response','Response','Response','Response','Resistance'),
				level = c('1','2','3A','4','R1'),
                 num = c(6,2,1,2,5))
re
re$Level <- factor(re$level, levels=c('1','2','3A','4','R1'), labels=c('1','2','3A','4','R1'))
re$drug <- factor(re$drug, levels=c('Response','Resistance'), labels=c('Response','Resistance'))
col = c(rgb(red=48,green=162,blue=72,alpha=255,max=255),
		rgb(red=22,green=121,blue=182,alpha=255,max=255),
		rgb(red=152,green=78,blue=158,alpha=255,max=255),
		rgb(red=67,green=67,blue=67,alpha=255,max=255),
		rgb(red=239,green=58,blue=38,alpha=255,max=255)
		)
p <- ggplot(data = re, aes(x = drug, y = num, fill = Level))+
     geom_bar(stat="identity",width=0.7) +
	 theme_classic() +
     scale_fill_manual(values = c('#30A248','#1679B6','#984E9E','#434343','#EF3A26')) +
	 xlab('') + ylab('Highest level of actionability\nin mCRC cases') +
	 theme(legend.position="bottom",
	 axis.text.y = element_text(color='black',size=14,margin=margin(0,5,0,0)),
	 legend.text = element_text(color='black',size=14,margin=margin(0,5,0,0)),
	 legend.title = element_text(color='black',size=14),
	 axis.title.y = element_text(color='black',size=14,margin=margin(0,10,0,0)),
	 axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=16,color='black'),
     axis.ticks = element_line(color='black'),
	 axis.line=element_line(color='black'),
     axis.title.x=element_text(color='black')) +
	 guides(fill=guide_legend(nrow=2, byrow=TRUE, title.position = "top"))
p
ggsave("OncoKB_actionability_barplot.pdf",p,width=2.7,height=5.5)
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
snv <- read.csv('NJMU_mCRC_OncoKB_mutation_forPlot.csv')
snv <- snv[,1:9]
snv$indsam <- gsub('_','-',paste(snv$Normal, snv$Class_a, sep=':'))
snv <- merge(snv, info[,c('indsam','sam_name')], by='indsam')
cnv <- read.csv('NJMU_mCRC_OncoKB_sCNA_forPlot.csv')
cnv <- subset(cnv, nchar(id)>0)
cnv <- cnv[,1:6]
sam <- read.csv('<PROJECT_ROOT>')
sam$id <- paste(sam$Tumor,sam$Normal,sep='_')
cnv <- merge(cnv, sam[,c('id','Normal','Class_a')], by='id')
cnv$indsam <- paste(cnv$Normal, cnv$Class_a, sep=':')
subset(cnv, !indsam %in% info$indsam)
cnv <- merge(cnv, info[,c('indsam','sam_name')], by='indsam')
snv <- snv[,c('Hugo_Symbol','Variant_Classification','level','sam_name')]
cnv <- cnv[,c('gene','Alterations','level','sam_name')]
intersect(snv$Hugo_Symbol, cnv$Hugo_Symbol)
dat <- rbind(snv,cnv)
dat$sam <- sapply(strsplit(dat$sam_name,' '),'[',1)
a <- subset(dat, level != 'R1')
a <- a[order(a$Hugo_Symbol, a$sam_name, a$Variant_Classification, a$level),]
a <- subset(a, !duplicated(paste(a$Hugo_Symbol, a$Variant_Classification, a$sam_name)))
r <- subset(dat, level == 'R1')
data <- rbind(a,r)
mut <- data
mut[grep("In_Frame",mut$Variant_Classification),'Variant_Classification'] = "In_frame_indel"
mut[grep("Frame_Shift",mut$Variant_Classification),'Variant_Classification'] = "Frameshift"
mut[grep("Missense_Mutation",mut$Variant_Classification),'Variant_Classification'] = "Missense"
mut[grep("Nonsense_Mutation",mut$Variant_Classification),'Variant_Classification'] = "Nonsense"
mut[grep("Splice_Site",mut$Variant_Classification),'Variant_Classification'] = "Splice_site"
mut[grep("Nonstop_Mutation",mut$Variant_Classification),'Variant_Classification'] = "Nonstop"
samid <- read.csv('../2_driverGene_waterfall_plot/plot/sample_id.csv')
gene_order <- c('KRAS','PIK3CA','SMARCB1','BRAF','CHEK1','EGFR','SMARCB1 del','FGFR1','MET','CCNE1','FGFR2','CDK4','MDM2')
maf_matrix <- matrix("" , ncol = length(unique(samid$id)) , nrow = length(unique(data$Hugo_Symbol)) , dimnames = list(gene_order, samid$id))
maf_matrix
for(gene in rownames(maf_matrix)){
		for(tumor in colnames(maf_matrix)){
			index <- which(mut$Hugo_Symbol==gene & mut$sam_name==tumor)
			if(length(index)==0){
				mutype=NA
			}else if(length(index)>=1){
				mutype <- paste0(unique( mut[index,'Variant_Classification'] ) ,collapse=';')
			}
			maf_matrix[ gene , tumor ] <-  mutype
			rm(mutype, index, tumor)
		}
		rm(gene)
	}
lev_ano <- matrix("" , ncol = length(unique(samid$id)) , nrow = length(unique(data$Hugo_Symbol)) , dimnames = list(gene_order, samid$id))
lev_ano
for(gene in rownames(lev_ano)){
		for(tumor in colnames(lev_ano)){
			index <- which(mut$Hugo_Symbol==gene & mut$sam_name==tumor)
			if(length(index)==0){
				mutype=NA
			}else if(length(index)>=1){
				mutype <- paste0(unique( mut[index,'level'] ) ,collapse=';\n')
			}
			lev_ano[ gene , tumor ] <-  mutype
			rm(mutype, index, tumor)
		}
		rm(gene)
	}
col = c(rgb(red=102,green=194,blue=164,alpha=255,max=255),
		rgb(red=231,green=137,blue=195,alpha=255,max=255),
		rgb(red=139,green=158,blue=201,alpha=255,max=255),
		rgb(red=234,green=173,blue=126,alpha=255,max=255),
		rgb(red=126,green=184,blue=215,alpha=255,max=255))
colors = structure(col, names = c('Missense','Frameshift','Splice_site','Amplification','Deletion'))
Heatmap(maf_matrix, name = "Mutation", col = colors, na_col = "#F2F2F2", rect_gp = gpar(col = "white", lwd = 2))
patient_id <- sapply(strsplit(colnames(maf_matrix),' '),'[',2)
Heatmap(maf_matrix, name = "Mutation", col = colors, na_col = "#F2F2F2",
	rect_gp = gpar(col = "white", lwd = 3),
	row_names_side = "left",
	column_names_side = "top",
	row_names_gp = gpar(fontface = "italic"),
	column_split = patient_id,
	column_title = NULL,
	column_gap = unit(2.5, "mm"),
	show_heatmap_legend = FALSE,
	cell_fun = function(j, i, x, y, width, height, fill) {
		if(!is.na(lev_ano[i, j])){
        grid.text(sprintf("%s", lev_ano[i, j]), x, y, gp = gpar(fontsize = 8))
		}
	}
		)
rownames(maf_matrix)[which(rownames(maf_matrix)=='SMARCB1 del ')] <- 'SMARCB1'
pdf("drug_target_heatmap.pdf",width=12, height=5.5)
Heatmap(maf_matrix, name = "Mutation", col = colors, na_col = "#F2F2F2",
	rect_gp = gpar(col = "white", lwd = 3),
	row_names_side = "left",
	column_names_side = "top",
	row_names_gp = gpar(fontface = "italic"),
	column_split = as.numeric(gsub('CRC','',patient_id)),
	row_split = c(1,1,1,1,1,2,2,2,2,2,2,2,2),
	column_title = NULL,
	row_title = NULL,
	column_gap = unit(2.5, "mm"),
	row_gap = unit(2.5, "mm"),
	show_heatmap_legend = FALSE,
	cell_fun = function(j, i, x, y, width, height, fill) {
		if(!is.na(lev_ano[i, j])){
        grid.text(sprintf("%s", lev_ano[i, j]), x, y, gp = gpar(fontsize = 8))
		}
	}
		)
dev.off()
ldg_hp = Legend(labels = names(colors) , title = "Mutation" ,
	 legend_gp = gpar(fill = colors , bar_width = 1 , fontsize = 12) , nrow = 1 ,
	 gap = unit(0.6, "cm"),
	 row_gap = unit(2, "mm")
	)
pdf('legend.pdf',height=2,width=12)
draw(ldg_hp,, x = unit(0.2, "npc"), y = unit(0.6, "npc"),just = c("left", "bottom"))
dev.off()
