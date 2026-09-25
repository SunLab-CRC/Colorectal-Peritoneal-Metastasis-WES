options(stringsAsFactors=F)
library(data.table)
setwd('<PROJECT_ROOT>')
rm(list=ls())
filename <- c(
'crc_msk_2017/MSK2017_mCRC.PERITONEUM.OMENTUM.ABDOMEN_primaries_Gene_Level_No.Mutation.csv',
'crc_msk_2017/MSK2017_mCRC.PERITONEUM.OMENTUM.ABDOMEN_Metastasis_Gene_Level_No.Mutation.csv',
'crc_msk_2017/MSK2017_Early_stage_primaries_Gene_Level_No.Mutation.csv',
'coadread_tcga_pan_can_atlas_2018/TCGA_COADREAD_StageI-III_MSS_Gene_Level_No.Mutation.csv',
'results_mutDriver_paired/mCRC_Gene_Level_No.Mutation.csv'
)
p1 <- read.csv(filename[5])
p2 <- read.csv(filename[1])
p3 <- read.csv(filename[2])
p <- merge(p1,p2,by='gene',all.y=T)
p <- merge(p,p3,by='gene')
p$mCRC_Met_T.njmu <- ifelse(is.na(p$mCRC_Met_T.njmu),0,p$mCRC_Met_T.njmu)
p$mCRC_Met_F.njmu <- ifelse(is.na(p$mCRC_Met_F.njmu),11,p$mCRC_Met_F.njmu)
m1 <- read.csv(filename[3])
m2 <- read.csv(filename[4])
dat <- merge(p,m1,by='gene')
dat <- merge(dat,m2,by='gene')
dat$mCRC_T <- dat$mCRC_Met_T.njmu + dat$mCRC_Met_T.msk + dat$mCRC_Met_T.mskpri
dat$mCRC_F <- dat$mCRC_Met_F.njmu + dat$mCRC_Met_F.msk + dat$mCRC_Met_F.mskpri
dat$Early_Stage_Primary_T <- dat$Early_Stage_Primary_T.msk + dat$Early_Stage_Primary_T.tcga
dat$Early_Stage_Primary_F <- dat$Early_Stage_Primary_F.msk + dat$Early_Stage_Primary_F.tcga
dat$pvalue <- NA
dat$OR <- NA
dat$L95 <- NA
dat$U95 <- NA
for(i in 1:nrow(dat)){
  a <- dat$Early_Stage_Primary_T[i]
  b <- dat$Early_Stage_Primary_F[i]
  d <- dat$mCRC_F[i]
  cc  <- dat$mCRC_T[i]
  mat <- matrix(c(cc,a,d,b),ncol=2,byrow=T)
  test <- fisher.test(mat,alternative = "greater")
  dat$pvalue[i] = test$p.value
  dat$OR[i] = test$estimate
  dat$L95[i] = test$conf.int[1]
  dat$U95[i] = test$conf.int[2]
  rm(a,b,cc,d,test,mat,i)
}
subset(dat, pvalue<0.05)
write.csv(dat,file='MSK_TCGA_results/Abdomen.NJMU.MSK.mCRC_vs_MSK.TCGA.Early.pri_Fisher.Test.csv',row.names=F,quote=F)
