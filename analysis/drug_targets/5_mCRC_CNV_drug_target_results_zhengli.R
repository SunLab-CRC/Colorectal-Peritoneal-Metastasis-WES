options(stringsAsFactors=F)
setwd('<PROJECT_ROOT>')
rm(list=ls())
dat <- read.csv('NJMU_mCRC_OncoKB_result.csv')
id <- read.csv('../0_patient_sample_info/SampleInformation.final.csv')
dat1 <- merge(dat[,7:ncol(dat)], id[,c('Tumor','Class_a')],by='Tumor')
data <- dat1
data$ps <- paste(data$Normal,data$Class_a)
snv <- subset(data, Variant_Classification == 'SNV')
re <- as.matrix(table(snv$ps, snv$Hugo_Symbol))
for(i in 1:nrow(snv)){
	re[snv$ps[i],snv$Hugo_Symbol[i]] <- paste(re[snv$ps[i],snv$Hugo_Symbol[i]], snv$Alterations[i], snv$Level[i], sep=';')
}
write.csv(re, file='OncoKB_SNV_result_table.csv',quote=F)
rm(re)
cnv <- subset(data, Variant_Classification == 'CNV')
re <- as.matrix(table(cnv$ps, cnv$Hugo_Symbol))
for(i in 1:nrow(cnv)){
	re[cnv$ps[i],cnv$Hugo_Symbol[i]] <- paste(re[cnv$ps[i],cnv$Hugo_Symbol[i]], cnv$Alterations[i], cnv$Level[i], sep=';')
}
write.csv(re, file='OncoKB_CNV_result_table.csv',quote=F)
rm(re)
snv <- read.csv('OncoKB_SNV_result_table.csv')
cnv <- read.csv('OncoKB_CNV_result_table.csv')
data <- merge(snv,cnv,by='X',all=T)
write.csv(data, file='OncoKB_ALL_result_table.csv',quote=F)
options(stringsAsFactors=F)
setwd('<PROJECT_ROOT>')
rm(list=ls())
dat <- read.csv('NJMU_mCRC_COSMIC_result.csv')
id <- read.csv('../0_patient_sample_info/SampleInformation.final.csv')
dat1 <- merge(dat, id[,c('Tumor','Class_a')],by='Tumor')
data <- dat1
data$ps <- paste(data$Normal,data$Class_a)
snv <- subset(data, Variant_Classification == 'SNV')
re <- as.matrix(table(snv$ps, snv$gene))
for(i in 1:nrow(snv)){
	re[snv$ps[i],snv$gene[i]] <- paste(re[snv$ps[i],snv$gene[i]], snv$mutation[i], snv$ACTIONABILITY_RANK[i], sep=';')
}
write.csv(re, file='COSMIC_SNV_result_table.csv',quote=F)
rm(re)
cnv <- subset(data, Variant_Classification == 'CNV')
re <- as.matrix(table(cnv$ps, cnv$gene))
for(i in 1:nrow(cnv)){
	re[cnv$ps[i],cnv$gene[i]] <- paste(re[cnv$ps[i],cnv$gene[i]], cnv$mutation[i], cnv$ACTIONABILITY_RANK[i], sep=';')
}
write.csv(re, file='COSMIC_CNV_result_table.csv',quote=F)
rm(re)
snv <- read.csv('COSMIC_SNV_result_table.csv')
cnv <- read.csv('COSMIC_CNV_result_table.csv')
data <- merge(snv,cnv,by='X',all=T)
write.csv(data, file='COSMIC_ALL_result_table.csv',quote=F)
