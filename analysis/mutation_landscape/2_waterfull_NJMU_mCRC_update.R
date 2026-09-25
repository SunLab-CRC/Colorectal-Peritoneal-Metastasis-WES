setwd('<PROJECT_ROOT>')
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(circlize)
rm(list=ls())
smg <- read.csv('../Colorectal_Cancer_Mut_Driver_Genes.csv')
info <- read.csv('<PROJECT_ROOT>')
info$indsam <- paste(info$Normal, info$Class_a, sep=':')
maf <- fread('NJMU_mCRC_ALL.Gene_mutation.maf',data.table=F)
maf$indsam <- paste(maf$Normal, maf$Class_a, sep=':')
dat <- merge(maf, info,by='indsam')
dat <- subset(dat, variant_allele_frequency>0)
count_mut <- function(x){ return ( length(unique(x)) ) }
mutn <- as.matrix (by(dat$mutation_id.x, dat$sam_name, count_mut) )
mutn  <- as.data.frame(mutn )
colnames(mutn ) <- 'n_mut'
mutn$sam_name <- rownames(mutn)
exon <- subset(dat, Variant_Classification %in% c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Silent','Splice_Site'))
count_mut <- function(x){ return ( length(unique(x)) ) }
mutn <- as.matrix (by(exon$mutation_id.x, exon$sam_name, count_mut) )
mutn  <- as.data.frame(mutn )
colnames(mutn ) <- 'n_mut'
mutn$sam_name <- rownames(mutn)
mut <- subset(dat , Variant_Classification %in% c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site') )
mut[grep("In_Frame",mut$Variant_Classification),'Variant_Classification'] = "In_frame_indel"
mut[grep("Frame_Shift",mut$Variant_Classification),'Variant_Classification'] = "Frameshift"
mut[grep("Missense_Mutation",mut$Variant_Classification),'Variant_Classification'] = "Missense"
mut[grep("Nonsense_Mutation",mut$Variant_Classification),'Variant_Classification'] = "Nonsense"
mut[grep("Splice_Site",mut$Variant_Classification),'Variant_Classification'] = "Splice_site"
mut[grep("Nonstop_Mutation",mut$Variant_Classification),'Variant_Classification'] = "Nonstop"
count_sample <- function(x){ return ( length(unique(x)) ) }
re <- as.matrix (by(mut$Normal.x, mut$Hugo_Symbol, count_sample) )
re <- as.data.frame(re)
colnames(re) <- 'n_sample'
re$gene <- rownames(re)
gene <- subset(re, (gene %in% smg$SYMBOL & n_sample>1) | n_sample > 2)
gene <- gene[order(gene$n_sample, decreasing=T),]
gene$freq <- gene$n_sample*100/11
gene_freq <- gene
dat <- subset(mut, Hugo_Symbol %in% gene$gene)
samid <- data.frame(id = unique(dat$sam_name))
samid$ind <- as.numeric( gsub('CRC','',sapply(strsplit(samid$id,' '),'[',1)) )
samid$sample_id <- sapply(strsplit(samid$id,' '),'[',2)
sample_id <- unique(samid$sample_id)
samid$sample_id <- factor(samid$sample_id, levels = c("P", "P1", "P2", "P3", "LN", "LN1", "LN2", "A", "A1", "A2", "A3", "A4", "A5", "A6"))
samid <- samid[order(samid$ind, samid$sample_id),]
write.csv(samid,file='sample_id.csv',row.names=F,quote=F)
maf_matrix <- matrix("" , ncol = length(unique(dat$sam_name)) , nrow = length(unique(dat$Hugo_Symbol)) , dimnames = list(gene_freq$gene, samid$id))
maf_matrix
for(gene in rownames(maf_matrix)){
		for(tumor in colnames(maf_matrix)){
			index <- which(dat$Hugo_Symbol==gene & dat$sam_name==tumor)
			if(length(index)==0){
				mutype=NA
			}else if(length(index)==1){
				mutype <- paste0(unique( dat[index,'Variant_Classification'] ) ,collapse=';')
			}else if(length(index)>1){
				mutype <- "Multiple_Hits"
			}
			maf_matrix[ gene , tumor ] <-  mutype
			rm(mutype, index, tumor)
		}
		rm(gene)
	}
driano <- rownames(maf_matrix) %in% smg$SYMBOL
driano[driano==TRUE] <- '*'
driano[driano==FALSE] <- ''
col = c(rgb(red=50,green=123,blue=186,alpha=255,max=255),
		rgb(red=222,green=10,blue=21,alpha=255,max=255),
		rgb(red=241,green=125,blue=0,alpha=255,max=255),
		rgb(red=179,green=213,blue=230,alpha=255,max=255),
		rgb(red=69,green=47,blue=140,alpha=255,max=255),
		rgb(red=147,green=12,blue=19,alpha=255,max=255),
		rgb(red=65,green=174,blue=119,alpha=255,max=255))
colors = structure(col, names = c('Missense','Nonsense','Frameshift','In_frame_indel','Splice_site','Nonstop','Multiple_Hits'))
Heatmap(maf_matrix, name = "Mutation", col = colors, na_col = "#F2F2F2", rect_gp = gpar(col = "white", lwd = 2))
patient_id <- sapply(strsplit(colnames(maf_matrix),' '),'[',1)
mutn <- mutn[match(samid$id,mutn$sam_name), ]
ht_opt("heatmap_row_names_gp" = gpar(fontface = "italic"),
    heatmap_border = FALSE,
    annotation_border = FALSE,
	COLUMN_ANNO_PADDING = unit(0.5, "cm"),
	ROW_ANNO_PADDING = unit(0.5, "cm")
)
col_ano <- default_axis_param("column")
col_ano$gp$fontsize <- 12
top_ano <- HeatmapAnnotation(
        No.mutation = anno_barplot(mutn$n_mut, height = unit(2.5, "cm"), bar_width = 0.7, gp = gpar(fill = "#B9BEC2", col = "#B9BEC2"), border=FALSE,axis_param=list(col_ano,
		gp=gpar(fontsize=12),
		at = c(0,50,100,150,200,250),
		labels = c('0','','100','','200','')))
		, annotation_name_rot = 90, annotation_name_side = "left" , border = FALSE)
row_ano <- default_axis_param("row")
row_ano$side <- "top"
row_ano$gp$fontsize <- 12
right_ano <- rowAnnotation(
        Number.of.affected.patients = anno_barplot(gene_freq$n_sample, width = unit(4, "cm"), bar_width = 0.7, gp = gpar(fill = "#B2D6E5", col = "#B2D6E5"), border=FALSE, annotation_name = "No. affected patients", axis_param=row_ano),annotation_name_rot = 0, annotation_name_side = "top")
Heatmap(maf_matrix, name = "Mutation", col = colors, na_col = "#F2F2F2",
	rect_gp = gpar(col = "white", lwd = 3),
	row_names_side = "left",
	row_names_gp = gpar(fontface = "italic"),
	top_annotation = top_ano,
	right_annotation = right_ano,
	column_split = patient_id,
	column_title = NULL,
	column_gap = unit(2.5, "mm"),
	show_heatmap_legend = FALSE
		)
pdf("Mut_WaterFall.new.pdf",width=11.5, height=7)
Heatmap(maf_matrix, name = "Mutation", col = colors, na_col = "#F2F2F2",
	rect_gp = gpar(col = "white", lwd = 3),
	row_names_side = "left",
	row_names_gp = gpar(fontface = "italic"),
	top_annotation = top_ano,
	column_split = as.numeric(gsub('CRC','',patient_id)),
	column_title = NULL,
	column_gap = unit(2.5, "mm"),
	show_heatmap_legend = FALSE
		)
dev.off()
