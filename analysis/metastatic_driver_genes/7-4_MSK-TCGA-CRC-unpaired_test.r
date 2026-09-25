cd <PROJECT_ROOT>
cp <PROJECT_ROOT> ./
unzip MSK_TCGA_20230817.zip
R
options(stringsAsFactors=F)
library(data.table)
setwd('<PROJECT_ROOT>')
rm(list=ls())
met <- read.delim("20230817/MSK_mCRC_Abdomen/scores.gistic")
met$chr <- paste('chr', met$Chromosome, sep='')
met <- met[,c('chr','Start','End','X.log10.q.value.','G.score','Type')]
write.table(met, file='MSK_Abdomen_Metastasis_seg.bed',row.names=F,col.names=F,quote=F,sep='\t')
rm(list=ls())
met <- read.delim("20230817/MSK_mCRC_Early/scores.gistic")
met$chr <- paste('chr', met$Chromosome, sep='')
met <- met[,c('chr','Start','End','X.log10.q.value.','G.score','Type')]
write.table(met, file='MSK_Early_Primary_seg.bed',row.names=F,col.names=F,quote=F,sep='\t')
rm(list=ls())
met <- read.delim("20230817/TCGA/scores.gistic")
met$chr <- paste('chr', met$Chromosome, sep='')
met <- met[,c('chr','Start','End','X.log10.q.value.','G.score','Type')]
write.table(met, file='MSK_TCGA_validation_Primary_seg.bed',row.names=F,col.names=F,quote=F,sep='\t')
cd <PROJECT_ROOT>
<PROJECT_ROOT> intersect \
-a <PROJECT_ROOT> \
-b MSK_Abdomen_Metastasis_seg.bed \
-wa -wb \
> MSK_Abdomen_Metastasis_gene_score.bed
wc -l MSK_Abdomen_Metastasis_gene_score.bed
<PROJECT_ROOT> intersect \
-a <PROJECT_ROOT> \
-b MSK_Early_Primary_seg.bed \
-wa -wb \
> MSK_Early_Primary_gene_score.bed
wc -l MSK_Early_Primary_gene_score.bed
<PROJECT_ROOT> intersect \
-a <PROJECT_ROOT> \
-b MSK_TCGA_validation_Primary_seg.bed \
-wa -wb \
> MSK_TCGA_Early_Primary_gene_score.bed
wc -l MSK_TCGA_Early_Primary_gene_score.bed
R
options(stringsAsFactors=F)
library(data.table)
library(MADAM)
setwd('<PROJECT_ROOT>')
rm(list=ls())
gene <- read.delim('<PROJECT_ROOT>',h=F)
gene <- subset(gene, V8 == 'AMP')
met <- read.delim("MSK_Abdomen_Metastasis_gene_score.bed",h=F)
met <- subset(met, V8 == 'AMP' & V14 == 'Amp')
umap <- subset(gene, ! V5 %in% met$V5)
umap$log10Q <- NA
umap$G.score <- NA
umap$Type <- NA
umap$chr.amp <- NA
umap$start.amp <- NA
umap$end.amp <- NA
seg <- read.table('MSK_Abdomen_Metastasis_seg.bed')
seg <- subset(seg, V6 == 'Amp')
for(i in 1:nrow(umap)){
	segs_u <- subset(seg, V1==umap$V1[i] & V2 > umap$V3[i])
	segs_u <- segs_u[which(segs_u$V2 == min(segs_u$V2)),]
	segs_l <- subset(seg, V1==umap$V1[i] & V3 < umap$V2[i])
	segs_l <- segs_l[which(segs_l$V3 == max(segs_l$V3)),]
	score <- mean(c(segs_l$V5, segs_u$V5))
	p <- fisher.method(matrix(10^(-1*c(segs_l$V4, segs_u$V4)), nrow=1))$p.value
	umap$log10Q[i] <- -1*log10(p)
	umap$G.score[i] <- score
	umap$Type[i] <- 'Amp'
	umap$start.amp[i] <- mean(segs_l$V2)
	umap$end.amp[i] <- mean(segs_u$V3)
	rm(segs_u, segs_l, score, p, i)
}
umap <- umap[,c('V1','V5','V6','log10Q','G.score','Type','start.amp','end.amp')]
dup <- subset(met,duplicated(met$V5))
dup <- subset(met, V5 %in% dup$V5)
re <- data.frame()
for(i in unique(dup$V5)){
	sdu <- subset(dup, V5 == i)
	score <- mean(sdu$V13)
	p <- fisher.method(matrix(10^(-1*sdu$V12), nrow=1))$p.value
	sre <- data.frame(ID=unique(sdu$V5),
	                  symbol = unique(sdu$V6),
					  log10Q = -1*log10(p),
					  G.score = score,
					  Type = 'Amp',
					  chr.amp = unique(sdu$V9),
					  start.amp = min(sdu$V10),
					  end.amp = max(sdu$V11)
					  )
	re <- rbind(re,sre)
	rm(sdu, score, p, i, sre)
}
ind <- subset(met, !V5 %in% dup$V5)
ind <- ind[,c('V5','V6','V12','V13','V14','V9','V10','V11')]
met_am <- rbind(umap, re, ind)
write.csv(met_am, file = 'MSK_Abdomen_Metastasis_SCNA.Amplification.Driver_Gscore.csv', row.names=F,quote=F)
system('head MSK_Abdomen_Metastasis_SCNA.Amplification.Driver_Gscore.csv')
system('wc -l MSK_Abdomen_Metastasis_SCNA.Amplification.Driver_Gscore.csv')
R
options(stringsAsFactors=F)
library(data.table)
library(MADAM)
setwd('<PROJECT_ROOT>')
rm(list=ls())
gene <- read.delim('<PROJECT_ROOT>',h=F)
gene <- subset(gene, V8 == 'AMP')
met <- read.delim("MSK_Early_Primary_gene_score.bed",h=F)
met <- subset(met, V8 == 'AMP' & V14 == 'Amp')
umap <- subset(gene, ! V5 %in% met$V5)
umap$log10Q <- NA
umap$G.score <- NA
umap$Type <- NA
umap$chr.amp <- NA
umap$start.amp <- NA
umap$end.amp <- NA
seg <- read.table('MSK_Early_Primary_seg.bed')
seg <- subset(seg, V6 == 'Amp')
for(i in 1:nrow(umap)){
	segs_u <- subset(seg, V1==umap$V1[i] & V2 > umap$V3[i])
	segs_u <- segs_u[which(segs_u$V2 == min(segs_u$V2)),]
	segs_l <- subset(seg, V1==umap$V1[i] & V3 < umap$V2[i])
	segs_l <- segs_l[which(segs_l$V3 == max(segs_l$V3)),]
	score <- mean(c(segs_l$V5, segs_u$V5))
	p <- fisher.method(matrix(10^(-1*c(segs_l$V4, segs_u$V4)), nrow=1))$p.value
	umap$log10Q[i] <- -1*log10(p)
	umap$G.score[i] <- score
	umap$Type[i] <- 'Amp'
	umap$start.amp[i] <- mean(segs_l$V2)
	umap$end.amp[i] <- mean(segs_u$V3)
	rm(segs_u, segs_l, score, p, i)
}
umap <- umap[,c('V1','V5','V6','log10Q','G.score','Type','start.amp','end.amp')]
dup <- subset(met,duplicated(met$V5))
dup <- subset(met, V5 %in% dup$V5)
re <- data.frame()
for(i in unique(dup$V5)){
	sdu <- subset(dup, V5 == i)
	score <- mean(sdu$V13)
	p <- fisher.method(matrix(10^(-1*sdu$V12), nrow=1))$p.value
	sre <- data.frame(ID=unique(sdu$V5),
	                  symbol = unique(sdu$V6),
					  log10Q = -1*log10(p),
					  G.score = score,
					  Type = 'Amp',
					  chr.amp = unique(sdu$V9),
					  start.amp = min(sdu$V10),
					  end.amp = max(sdu$V11)
					  )
	re <- rbind(re,sre)
	rm(sdu, score, p, i, sre)
}
ind <- subset(met, !V5 %in% dup$V5)
ind <- ind[,c('V5','V6','V12','V13','V14','V9','V10','V11')]
met_am <- rbind(umap, re, ind)
write.csv(met_am, file = 'MSK_Early_Primary_SCNA.Amplification.Driver_Gscore.csv', row.names=F,quote=F)
system('head MSK_Early_Primary_SCNA.Amplification.Driver_Gscore.csv')
system('wc -l MSK_Early_Primary_SCNA.Amplification.Driver_Gscore.csv')
R
options(stringsAsFactors=F)
library(data.table)
library(MADAM)
setwd('<PROJECT_ROOT>')
rm(list=ls())
gene <- read.delim('<PROJECT_ROOT>',h=F)
gene <- subset(gene, V8 == 'AMP')
met <- read.delim("MSK_TCGA_Early_Primary_gene_score.bed",h=F)
met <- subset(met, V8 == 'AMP' & V14 == 'Amp')
umap <- subset(gene, ! V5 %in% met$V5)
umap$log10Q <- NA
umap$G.score <- NA
umap$Type <- NA
umap$chr.amp <- NA
umap$start.amp <- NA
umap$end.amp <- NA
seg <- read.table('MSK_TCGA_validation_Primary_seg.bed')
seg <- subset(seg, V6 == 'Amp')
for(i in 1:nrow(umap)){
	segs_u <- subset(seg, V1==umap$V1[i] & V2 > umap$V3[i])
	segs_u <- segs_u[which(segs_u$V2 == min(segs_u$V2)),]
	segs_l <- subset(seg, V1==umap$V1[i] & V3 < umap$V2[i])
	segs_l <- segs_l[which(segs_l$V3 == max(segs_l$V3)),]
	score <- mean(c(segs_l$V5, segs_u$V5))
	p <- fisher.method(matrix(10^(-1*c(segs_l$V4, segs_u$V4)), nrow=1))$p.value
	umap$log10Q[i] <- -1*log10(p)
	umap$G.score[i] <- score
	umap$Type[i] <- 'Amp'
	umap$start.amp[i] <- mean(segs_l$V2)
	umap$end.amp[i] <- mean(segs_u$V3)
	rm(segs_u, segs_l, score, p, i)
}
umap <- umap[,c('V1','V5','V6','log10Q','G.score','Type','start.amp','end.amp')]
dup <- subset(met,duplicated(met$V5))
dup <- subset(met, V5 %in% dup$V5)
re <- data.frame()
for(i in unique(dup$V5)){
	sdu <- subset(dup, V5 == i)
	score <- mean(sdu$V13)
	p <- fisher.method(matrix(10^(-1*sdu$V12), nrow=1))$p.value
	sre <- data.frame(ID=unique(sdu$V5),
	                  symbol = unique(sdu$V6),
					  log10Q = -1*log10(p),
					  G.score = score,
					  Type = 'Amp',
					  chr.amp = unique(sdu$V9),
					  start.amp = min(sdu$V10),
					  end.amp = max(sdu$V11)
					  )
	re <- rbind(re,sre)
	rm(sdu, score, p, i, sre)
}
ind <- subset(met, !V5 %in% dup$V5)
ind <- ind[,c('V5','V6','V12','V13','V14','V9','V10','V11')]
met_am <- rbind(umap, re, ind)
write.csv(met_am, file = 'MSK_TCGA_Early_Primary_SCNA.Amplification.Driver_Gscore.csv', row.names=F,quote=F)
system('head MSK_TCGA_Early_Primary_SCNA.Amplification.Driver_Gscore.csv')
system('wc -l MSK_TCGA_Early_Primary_SCNA.Amplification.Driver_Gscore.csv')
R
options(stringsAsFactors=F)
library(data.table)
library(MADAM)
setwd('<PROJECT_ROOT>')
rm(list=ls())
gene <- read.delim('<PROJECT_ROOT>',h=F)
gene <- subset(gene, V8 == 'DEL')
met <- read.delim("MSK_Abdomen_Metastasis_gene_score.bed",h=F)
met <- subset(met, V8 == 'DEL' & V14 == 'Del')
umap <- subset(gene, ! V5 %in% met$V5)
umap$log10Q <- NA
umap$G.score <- NA
umap$Type <- NA
umap$chr.del <- NA
umap$start.del <- NA
umap$end.del <- NA
seg <- read.table('MSK_Abdomen_Metastasis_seg.bed')
seg <- subset(seg, V6 == 'Del')
for(i in 1:nrow(umap)){
	segs_u <- subset(seg, V1==umap$V1[i] & V2 > umap$V3[i])
	segs_u <- segs_u[which(segs_u$V2 == min(segs_u$V2)),]
	segs_l <- subset(seg, V1==umap$V1[i] & V3 < umap$V2[i])
	segs_l <- segs_l[which(segs_l$V3 == max(segs_l$V3)),]
	score <- mean(c(segs_l$V5, segs_u$V5))
	p <- fisher.method(matrix(10^(-1*c(segs_l$V4, segs_u$V4)), nrow=1))$p.value
	umap$log10Q[i] <- -1*log10(p)
	umap$G.score[i] <- score
	umap$Type[i] <- 'Del'
	umap$start.del[i] <- mean(segs_l$V2)
	umap$end.del[i] <- mean(segs_u$V3)
	rm(segs_u, segs_l, score, p, i)
}
umap <- umap[,c('V1','V5','V6','log10Q','G.score','Type','start.del','end.del')]
dup <- subset(met,duplicated(met$V5))
dup <- subset(met, V5 %in% dup$V5)
re <- data.frame()
for(i in unique(dup$V5)){
	sdu <- subset(dup, V5 == i)
	score <- mean(sdu$V13)
	p <- fisher.method(matrix(10^(-1*sdu$V12), nrow=1))$p.value
	sre <- data.frame(ID=unique(sdu$V5),
	                  symbol = unique(sdu$V6),
					  log10Q = -1*log10(p),
					  G.score = score,
					  Type = 'Del',
					  chr.del = unique(sdu$V9),
					  start.del = min(sdu$V10),
					  end.del = max(sdu$V11)
					  )
	re <- rbind(re,sre)
	rm(sdu, score, p, i, sre)
}
ind <- subset(met, !V5 %in% dup$V5)
ind <- ind[,c('V5','V6','V12','V13','V14','V9','V10','V11')]
met_am <- rbind(umap, re, ind)
write.csv(met_am, file = 'MSK_Abdomen_Metastasis_SCNA.Deletion.Driver_Gscore.csv', row.names=F,quote=F)
system('head MSK_Abdomen_Metastasis_SCNA.Deletion.Driver_Gscore.csv')
system('wc -l MSK_Abdomen_Metastasis_SCNA.Deletion.Driver_Gscore.csv')
R
options(stringsAsFactors=F)
library(data.table)
library(MADAM)
setwd('<PROJECT_ROOT>')
rm(list=ls())
gene <- read.delim('<PROJECT_ROOT>',h=F)
gene <- subset(gene, V8 == 'DEL')
met <- read.delim("MSK_Early_Primary_gene_score.bed",h=F)
met <- subset(met, V8 == 'DEL' & V14 == 'Del')
umap <- subset(gene, ! V5 %in% met$V5)
umap$log10Q <- NA
umap$G.score <- NA
umap$Type <- NA
umap$chr.del <- NA
umap$start.del <- NA
umap$end.del <- NA
seg <- read.table('MSK_Early_Primary_seg.bed')
seg <- subset(seg, V6 == 'Del')
for(i in 1:nrow(umap)){
	segs_u <- subset(seg, V1==umap$V1[i] & V2 > umap$V3[i])
	segs_u <- segs_u[which(segs_u$V2 == min(segs_u$V2)),]
	segs_l <- subset(seg, V1==umap$V1[i] & V3 < umap$V2[i])
	segs_l <- segs_l[which(segs_l$V3 == max(segs_l$V3)),]
	score <- mean(c(segs_l$V5, segs_u$V5))
	p <- fisher.method(matrix(10^(-1*c(segs_l$V4, segs_u$V4)), nrow=1))$p.value
	umap$log10Q[i] <- -1*log10(p)
	umap$G.score[i] <- score
	umap$Type[i] <- 'Del'
	umap$start.del[i] <- mean(segs_l$V2)
	umap$end.del[i] <- mean(segs_u$V3)
	rm(segs_u, segs_l, score, p, i)
}
umap <- umap[,c('V1','V5','V6','log10Q','G.score','Type','start.del','end.del')]
dup <- subset(met,duplicated(met$V5))
dup <- subset(met, V5 %in% dup$V5)
re <- data.frame()
for(i in unique(dup$V5)){
	sdu <- subset(dup, V5 == i)
	score <- mean(sdu$V13)
	p <- fisher.method(matrix(10^(-1*sdu$V12), nrow=1))$p.value
	sre <- data.frame(ID=unique(sdu$V5),
	                  symbol = unique(sdu$V6),
					  log10Q = -1*log10(p),
					  G.score = score,
					  Type = 'Del',
					  chr.del = unique(sdu$V9),
					  start.del = min(sdu$V10),
					  end.del = max(sdu$V11)
					  )
	re <- rbind(re,sre)
	rm(sdu, score, p, i, sre)
}
ind <- subset(met, !V5 %in% dup$V5)
ind <- ind[,c('V5','V6','V12','V13','V14','V9','V10','V11')]
met_am <- rbind(umap, re, ind)
write.csv(met_am, file = 'MSK_Early_Primary_SCNA.Deletion.Driver_Gscore.csv', row.names=F,quote=F)
system('head MSK_Early_Primary_SCNA.Deletion.Driver_Gscore.csv')
system('wc -l MSK_Early_Primary_SCNA.Deletion.Driver_Gscore.csv')
R
options(stringsAsFactors=F)
library(data.table)
library(MADAM)
setwd('<PROJECT_ROOT>')
rm(list=ls())
gene <- read.delim('<PROJECT_ROOT>',h=F)
gene <- subset(gene, V8 == 'DEL')
met <- read.delim("MSK_TCGA_Early_Primary_gene_score.bed",h=F)
met <- subset(met, V8 == 'DEL' & V14 == 'Del')
umap <- subset(gene, ! V5 %in% met$V5)
umap$log10Q <- NA
umap$G.score <- NA
umap$Type <- NA
umap$chr.del <- NA
umap$start.del <- NA
umap$end.del <- NA
seg <- read.table('MSK_TCGA_validation_Primary_seg.bed')
seg <- subset(seg, V6 == 'Del')
for(i in 1:nrow(umap)){
	segs_u <- subset(seg, V1==umap$V1[i] & V2 > umap$V3[i])
	segs_u <- segs_u[which(segs_u$V2 == min(segs_u$V2)),]
	segs_l <- subset(seg, V1==umap$V1[i] & V3 < umap$V2[i])
	segs_l <- segs_l[which(segs_l$V3 == max(segs_l$V3)),]
	score <- mean(c(segs_l$V5, segs_u$V5))
	p <- fisher.method(matrix(10^(-1*c(segs_l$V4, segs_u$V4)), nrow=1))$p.value
	umap$log10Q[i] <- -1*log10(p)
	umap$G.score[i] <- score
	umap$Type[i] <- 'Del'
	umap$start.del[i] <- mean(segs_l$V2)
	umap$end.del[i] <- mean(segs_u$V3)
	rm(segs_u, segs_l, score, p, i)
}
umap <- umap[,c('V1','V5','V6','log10Q','G.score','Type','start.del','end.del')]
dup <- subset(met,duplicated(met$V5))
dup <- subset(met, V5 %in% dup$V5)
re <- data.frame()
for(i in unique(dup$V5)){
	sdu <- subset(dup, V5 == i)
	score <- mean(sdu$V13)
	p <- fisher.method(matrix(10^(-1*sdu$V12), nrow=1))$p.value
	sre <- data.frame(ID=unique(sdu$V5),
	                  symbol = unique(sdu$V6),
					  log10Q = -1*log10(p),
					  G.score = score,
					  Type = 'Del',
					  chr.del = unique(sdu$V9),
					  start.del = min(sdu$V10),
					  end.del = max(sdu$V11)
					  )
	re <- rbind(re,sre)
	rm(sdu, score, p, i, sre)
}
ind <- subset(met, !V5 %in% dup$V5)
ind <- ind[,c('V5','V6','V12','V13','V14','V9','V10','V11')]
met_am <- rbind(umap, re, ind)
write.csv(met_am, file = 'MSK_TCGA_Early_Primary_SCNA.Deletion.Driver_Gscore.csv', row.names=F,quote=F)
system('head MSK_TCGA_Early_Primary_SCNA.Deletion.Driver_Gscore.csv')
system('wc -l MSK_TCGA_Early_Primary_SCNA.Deletion.Driver_Gscore.csv')
cd <PROJECT_ROOT>
cp MSK_Abdomen_Metastasis_SCNA.Amplification.Driver_Gscore.csv <PROJECT_ROOT>
cp MSK_Early_Primary_SCNA.Amplification.Driver_Gscore.csv <PROJECT_ROOT>
cp MSK_TCGA_Early_Primary_SCNA.Amplification.Driver_Gscore.csv <PROJECT_ROOT>
cp MSK_Abdomen_Metastasis_SCNA.Deletion.Driver_Gscore.csv <PROJECT_ROOT>
cp MSK_Early_Primary_SCNA.Deletion.Driver_Gscore.csv <PROJECT_ROOT>
cp MSK_TCGA_Early_Primary_SCNA.Deletion.Driver_Gscore.csv <PROJECT_ROOT>
