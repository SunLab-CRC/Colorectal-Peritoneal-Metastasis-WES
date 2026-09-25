R
options(stringsAsFactors=F)
setwd('<PROJECT_ROOT>')
library(data.table)
rm(list=ls())
work_dir <- "<PROJECT_ROOT>"
seg_file <- "<PROJECT_ROOT>"
sample_file <- paste(work_dir,"/config/tumor_normal.list",sep="")
class_order <- paste(work_dir,"/config/Class_order.list",sep="")
class_a_order <- paste(work_dir,"/config/Class_a_order.list",sep="")
info<- data.frame(fread(sample_file))
class_order <- data.frame(fread(class_order))
class_a_order <- data.frame(fread(class_a_order))
seg <- data.frame(fread(seg_file))
info$ID <- paste0(info$Tumor , "_" , info$Normal )
dat <- merge(info,seg,by.x="ID" , by.y="Sample")
dat$Class <- factor(dat$Class,levels= unique(class_order$Class), ordered=TRUE)
dat$Class_a <- factor(dat$Class_a,levels= unique(class_a_order$Class), ordered=TRUE)
min_burden <- -0.01
max_burden <- 1
t_gain <- log2(2.5/2)
t_loss <- log2(1.5/2)
dat <- subset( dat , Chromosome!="X")
dat$TCN <- 2*2^dat$Corrected_logR
dat$Length <- as.numeric(abs(dat$Start - dat$End))
loss_region <- dat[which(dat$Corrected_logR <= t_loss),]
gain_region <- dat[which(dat$Corrected_logR >= t_gain),]
result_loss <- c()
result_gain <- c()
result_cna <- c()
for(Tumor in unique(dat$Tumor)){
	class_s <- unique(dat[which(dat$Tumor==Tumor),'Class_a'])
	normal <- unique(dat[which(dat$Tumor==Tumor),'Normal'])
	loss_rate <- sum(loss_region[which(loss_region$Tumor==Tumor),'Length'])/sum(dat[which(dat$Tumor==Tumor),'Length'])
	result_loss <- rbind(result_loss,data.frame(Tumor=Tumor,Normal=normal,Rate=loss_rate,Class=class_s))
	gain_rate <- sum(gain_region[which(gain_region$Tumor==Tumor),'Length'])/sum(dat[which(dat$Tumor==Tumor),'Length'])
	result_gain <- rbind(result_gain,data.frame(Tumor=Tumor,Normal=normal,Rate=gain_rate,Class=class_s))
	cna_rate <- (sum(loss_region[which(loss_region$Tumor==Tumor),'Length']) + sum(gain_region[which(gain_region$Tumor==Tumor),'Length']))/sum(dat[which(dat$Tumor==Tumor),'Length'])
	result_cna <- rbind(result_cna,data.frame(Tumor=Tumor,Normal=normal,Rate=cna_rate,Class=class_s))
}
dir.create("<PROJECT_ROOT>")
write.csv(result_cna, file='<PROJECT_ROOT>',row.names=F,quote=F)
write.csv(result_loss, file='<PROJECT_ROOT>',row.names=F,quote=F)
write.csv(result_gain, file='<PROJECT_ROOT>',row.names=F,quote=F)
options(stringsAsFactors=F)
setwd('<PROJECT_ROOT>')
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
rm(list=ls())
sam <- read.csv('sCNA_AmpDel_burden.csv')
id <- read.csv('../../0_patient_sample_info/Clonevol_sample_used.csv')
sam <- subset(sam, paste(Normal,Class) %in% paste(id$Normal,id$Class_a))
sam$Class_a <- sam$Class
sam$Class <- sapply(strsplit(sam$Class,'-'),'[',1)
dat <- data.frame()
for(i in unique(sam$Normal)){
	sams <- subset(sam, Normal == i)
	p <- subset(sams, Class == 'Primary')
	primary_mean <- mean(p$Rate)
	rm(p)
	n <- subset(sams, Class == 'Nodule')
	Nodule_mean <- mean(n$Rate)
	rm(n)
	l <- subset(sams, Class == 'Lymph')
	if(nrow(l)==0){
	  Lymph_mean <- NA
	}else{
	  Lymph_mean <- mean(l$Rate)
	}
	rm(l)
	m <- subset(sams, Class %in% c('Lymph','Nodule'))
	Met_mean <- mean(m$Rate)
	rm(m)
	dat <- rbind(dat,
	           data.frame(id = i,
			              primary_mean = primary_mean,
						  Nodule_mean  = Nodule_mean,
						  Lymph_mean   = Lymph_mean,
						  Met_mean     = Met_mean
						  ))
}
wilcox.test(dat$primary_mean, dat$Nodule_mean, paired = TRUE)
cor.test(dat$primary_mean, dat$Nodule_mean)
library(ggplot2)
library(ggprism)
sda <- rbind(data.frame(id=dat$id, rate=dat$primary_mean, Class='P'),
  data.frame(id=dat$id, rate=dat$Nodule_mean, Class='N'),
  data.frame(id=dat$id, rate=dat$Lymph_mean, Class='LN')
  )
sda <- subset(sda, !is.na(rate))
setheme <- theme(axis.text.y = element_text(size=16, color = "black"),
		  axis.text.x = element_text(size=16, color = "black"),
		  axis.title.y = element_text(size=16, color = "black", face = "bold"),
		  axis.title.x = element_text(size=16, color = "black", face = "bold"),
		  plot.title = element_text(size=16, color = "black",hjust = 0.5,margin = margin(t = 1, r = 0, b = 0.5, l = 0,"cm"), , face = "bold")) +
	theme(axis.line = element_line(colour = 'black', size = 0.5)) +
	theme(plot.margin = margin(t=0.2, r=2, b=0.2, l=1, "cm"))
p <- subset(sda,Class=='P')
n <- subset(sda,Class=='N')
pn <- merge(p,n,by='id')
p0 <- ggplot(aes(x = rate.x, y = rate.y),data=pn) +
      geom_point(fill='#3F76B4',color='#212121',alpha=0.6, size=4, shape=21) +
	  theme_classic() +
	  ylab("sCNA burden in\nabdominal metastases") +
	  xlab("sCNA burden in primary CRC") +
	  ggtitle("   ") +
	  stat_smooth(method = "lm",formula = y ~ x,colour = "lightgrey",fill='lightgrey',alpha=0.6) +
	  setheme
p1 <- p0 +
      geom_text(x = 0.2, y = 0.8,label=expression(paste(r^2," = 0.70",sep='')), size=5) +
	  geom_text(x = 0.2, y = 0.7,label=expression(paste(italic(P)," = 0.02",sep='')), size=5)
p1
ggsave('sCNA.rate_Primary_Nodule_scatterPlot.pdf',p1, width=5 ,height= 3.8)
setheme <- theme(axis.text.y = element_text(size=16, color = "black"),
		  axis.text.x = element_text(size=16, color = "black"),
		  axis.title.y = element_text(size=16, color = "black", face = "bold"),
		  axis.title.x = element_text(size=16, color = "black", face = "bold"),
		  plot.title = element_text(size=16, color = "black",hjust = 0.5,margin = margin(t = 1, r = 0, b = 0.5, l = 0,"cm"), , face = "bold")) +
	theme(axis.line = element_line(colour = 'black', size = 0.5)) +
	theme(plot.margin = margin(t=0.2, r=2, b=0.2, l=1, "cm"))
sda$Type <- ifelse(sda$Class=='P','P\n(n=11)',ifelse(sda$Class=='LN','LN\n(n=3)',ifelse(sda$Class=='N','A\n(n=11)',NA)))
sda$Type <- factor(sda$Type,levels=c('P\n(n=11)','A\n(n=11)','LN\n(n=3)'),labels=c('P\n(n=11)','A\n(n=11)','LN\n(n=3)'))
p0 <- ggplot(aes(x=Type,y=rate),data=sda) +
      geom_boxplot(fill='#3F76B4',color='#212121',width = 0.6,alpha=0.2,outlier.shape = NA) +
	  geom_jitter(width=0.1,color='#212121',size=3.5,alpha=0.3) +
	  theme_classic() +
	  xlab('') +
	  ylab('sCNA burden') +
	  ylim(0.1,1.2) +
	  setheme
wilcox.test(dat$primary_mean, dat$Nodule_mean, paired = TRUE)
wilcox.test(dat$primary_mean, dat$Lymph_mean, paired = TRUE)
p_vals <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~label,
  "P\n(n=11)",   "A\n(n=11)",     0.32, 0.9, "p.exprs.ital",
  "P\n(n=11)",   "LN\n(n=3)",     0.75, 1.1, "p.exprs.ital"
)
p1 = p0 + add_pvalue(p_vals, label = "P = {p.adj}", tip.length = 0.02, label.size = 5)
p1
ggsave('sCNA_Primary.VS.Met_paired.pdf',p1, width=5 ,height= 3.6)
