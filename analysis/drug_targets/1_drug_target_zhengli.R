setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
dat <- read.delim('OncoKB_Level1.txt',h=F)
re <- data.frame()
for(i in 1:(nrow(dat)/4)){
	sre <- data.frame(gene = dat$V1[4*(i-1)+1],
	                  alt = dat$V1[4*(i-1)+2],
					  cancerType = dat$V1[4*(i-1)+3],
					  drugs = dat$V1[4*(i-1)+4])
	re <- rbind(re,sre)
	rm(sre, i)
}
write.table(re,file='OncoKB_Level_1.txt',row.names=F,quote=F,sep='$')
rm(list=ls())
dat <- read.delim('OncoKB_Level2.txt',h=F)
re <- data.frame()
for(i in 1:(nrow(dat)/4)){
	sre <- data.frame(gene = dat$V1[4*(i-1)+1],
	                  alt = dat$V1[4*(i-1)+2],
					  cancerType = dat$V1[4*(i-1)+3],
					  drugs = dat$V1[4*(i-1)+4])
	re <- rbind(re,sre)
	rm(sre, i)
}
write.table(re,file='OncoKB_Level_2.txt',row.names=F,quote=F,sep='$')
rm(list=ls())
dat <- read.delim('OncoKB_Level3.txt',h=F)
re <- data.frame()
for(i in 1:(nrow(dat)/4)){
	sre <- data.frame(gene = dat$V1[4*(i-1)+1],
	                  alt = dat$V1[4*(i-1)+2],
					  cancerType = dat$V1[4*(i-1)+3],
					  drugs = dat$V1[4*(i-1)+4])
	re <- rbind(re,sre)
	rm(sre, i)
}
write.table(re,file='OncoKB_Level_3.txt',row.names=F,quote=F,sep='$')
rm(list=ls())
dat <- read.delim('OncoKB_Level4.txt',h=F)
re <- data.frame()
for(i in 1:(nrow(dat)/4)){
	sre <- data.frame(gene = dat$V1[4*(i-1)+1],
	                  alt = dat$V1[4*(i-1)+2],
					  cancerType = dat$V1[4*(i-1)+3],
					  drugs = dat$V1[4*(i-1)+4])
	re <- rbind(re,sre)
	rm(sre, i)
}
write.table(re,file='OncoKB_Level_4.txt',row.names=F,quote=F,sep='$')
rm(list=ls())
dat <- read.delim('OncoKB_R1.txt',h=F)
re <- data.frame()
for(i in 1:(nrow(dat)/4)){
	sre <- data.frame(gene = dat$V1[4*(i-1)+1],
	                  alt = dat$V1[4*(i-1)+2],
					  cancerType = dat$V1[4*(i-1)+3],
					  drugs = dat$V1[4*(i-1)+4])
	re <- rbind(re,sre)
	rm(sre, i)
}
write.table(re,file='OncoKB_R1.txt',row.names=F,quote=F,sep='$')
rm(list=ls())
dat <- read.delim('OncoKB_R2.txt',h=F)
re <- data.frame()
for(i in 1:(nrow(dat)/4)){
	sre <- data.frame(gene = dat$V1[4*(i-1)+1],
	                  alt = dat$V1[4*(i-1)+2],
					  cancerType = dat$V1[4*(i-1)+3],
					  drugs = dat$V1[4*(i-1)+4])
	re <- rbind(re,sre)
	rm(sre, i)
}
write.table(re,file='OncoKB_R2.txt',row.names=F,quote=F,sep='$')
setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
a1 <- read.table('OncoKB_Level_1.txt',sep='$',h=T)
a1$TxLevel <- 'Level_1'
a2 <- read.table('OncoKB_Level_2.txt',sep='$',h=T)
a2$TxLevel <- 'Level_2'
a3 <- read.table('OncoKB_Level_3.txt',sep='$',h=T)
a3$TxLevel <- 'Level_3'
a4 <- read.table('OncoKB_Level_4.txt',sep='$',h=T)
a4$TxLevel <- 'Level_4'
a5 <- read.table('OncoKB_R1.txt',sep='$',h=T)
a5$TxLevel <- 'Level_5'
a6 <- read.table('OncoKB_R2.txt',sep='$',h=T)
a6$TxLevel <- 'Level_6'
gen <- rbind(a1,a2,a3,a4,a5,a6)
write.table(gen,file='OncoKB_data_all.txt',row.names=F,quote=F,sep='$')
rm(list=ls())
dat <- read.delim('Actionability_AllData_v9_GRCh37.tsv')
dat <- subset(dat,!duplicated(dat[,1:5]))
write.csv(dat[,1:5],file='COSMIC_unique_gene.csv',row.names=F,quote=F)
re <- data.frame()
for(i in 1:(nrow(dat)/4)){
	sre <- data.frame(gene = dat$V1[4*(i-1)+1],
	                  alt = dat$V1[4*(i-1)+2],
					  cancerType = dat$V1[4*(i-1)+3],
					  drugs = dat$V1[4*(i-1)+4])
	re <- rbind(re,sre)
	rm(sre, i)
}
write.csv(re,file='OncoKB_R2.csv',row.names=F,quote=F)
