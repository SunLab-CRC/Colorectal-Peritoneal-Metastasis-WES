setwd('<PROJECT_ROOT>')
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
g <- read.delim('intOGen_release_date_2020.02.01/IntOGen-Cohorts-20200201/2020-02-02_IntOGen-Cohorts-20200213/cohorts.tsv')
g <- subset(g, !duplicated(paste(CANCER_TYPE, CANCER_TYPE_NAME)))
write.csv(g[,c('CANCER_TYPE','CANCER_TYPE_NAME')],file='IntOGen_Cancer_Type.csv',row.names=F,quote=F)
rm(g)
ing <- read.delim('intOGen_release_date_2020.02.01/IntOGen-Drivers-20200201/2020-02-02_IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv')
ing <- subset(ing, CANCER_TYPE %in% c("COREAD"))
ing <- data.frame(SYMBOL = unique(ing$SYMBOL), database = "IntOGen-COREAD")
tcga12 <- read.delim('2012_Nature_TCGA_CRC/TCGA_Nature_2012.txt')
tcga18 <- read.csv("2018_Cell_TCGA_PanCancer/TCGA_PanCancer_Cell2018.csv")
tcga18 <- subset(tcga18, Cancer == "COADREAD")
tcga18 <- tcga18[,1,drop=F]
ing$Cancer.Subtype <- "COREAD"
colnames(tcga18)[1] <- 'SYMBOL'
tcga18$Cancer.Subtype <- "COREAD"
tcga18$database <- "TCGA_Cell.2018"
dat <- rbind(ing, tcga12, tcga18)
dat$SYMBOL <- gsub(' ','',dat$SYMBOL)
cosm <- read.csv("COSMIC_CancerGeneCensus_v97/COSMIC_colorectal_cancer.csv",h=F)
nodri <- read.csv('COSMIC_CancerGeneCensus_v97/COSMIC_nonMutDriver.csv')
cosm <- subset(cosm, !V1 %in% nodri$GeneSymbol)
cosm <- cosm[,1,drop=F]
colnames(cosm)[1] <- 'SYMBOL'
cosm$Cancer.Subtype <- "colorectal cancer"
cosm$database <- "COSMIC"
data <- rbind(dat, cosm)
data$SYMBOL <- gsub(' ','',data$SYMBOL)
write.csv(data,file='Colorectal_Cancer_Mut_Driver_Genes.csv',row.names=F,quote=F)
