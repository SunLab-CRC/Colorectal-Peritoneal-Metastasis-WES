options(stringsAsFactors=F)
setwd('<PROJECT_ROOT>')
library(ggplot2)
library(ggprism)
rm(list=ls())
sam <- read.csv('SampleInformation.final.csv')
sam <- subset(sam, !is.na(Purity) & !is.na(Ploidy))
id <- read.csv('../Clonevol_sample_used.csv')
sam <- subset(sam, paste(Normal,Class_a) %in% paste(id$Normal,id$Class_a))
dat <- data.frame()
for(i in unique(sam$Normal)){
	sams <- subset(sam, Normal == i)
	p <- subset(sams, Class == 'Primary')
	Ploidy_primary_mean <- mean(p$Ploidy)
	Purity_primary_mean <- mean(p$Purity)
	rm(p)
	n <- subset(sams, Class == 'Nodule')
	Ploidy_Nodule_mean <- mean(n$Ploidy)
	Purity_Nodule_mean <- mean(n$Purity)
	rm(n)
	l <- subset(sams, Class == 'Lymph')
	if(nrow(l)==0){
	  Ploidy_Lymph_mean <- NA
	  Purity_Lymph_mean <- NA
	}else{
	  Ploidy_Lymph_mean <- mean(l$Ploidy)
	  Purity_Lymph_mean <- mean(l$Purity)
	}
	rm(l)
	m <- subset(sams, Class %in% c('Lymph','Nodule'))
	Ploidy_Met_mean <- mean(m$Ploidy)
	Purity_Met_mean <- mean(m$Purity)
	rm(m)
	dat <- rbind(dat,
	           data.frame(id = i,
			              Ploidy_primary_mean = Ploidy_primary_mean,
						  Ploidy_Nodule_mean  = Ploidy_Nodule_mean,
						  Ploidy_Lymph_mean   = Ploidy_Lymph_mean,
						  Ploidy_Met_mean     = Ploidy_Met_mean,
						  Purity_primary_mean = Purity_primary_mean,
						  Purity_Nodule_mean  = Purity_Nodule_mean,
						  Purity_Lymph_mean   = Purity_Lymph_mean,
						  Purity_Met_mean     = Purity_Met_mean
						  ))
}
wilcox.test(dat$Ploidy_primary_mean, dat$Ploidy_Nodule_mean, paired = TRUE)
cor.test(dat$Ploidy_primary_mean, dat$Ploidy_Nodule_mean, method = "spearman")
wilcox.test(dat$Ploidy_primary_mean, dat$Ploidy_Lymph_mean, paired = TRUE)
cor.test(dat$Ploidy_primary_mean, dat$Ploidy_Lymph_mean, method = "spearman")
wilcox.test(dat$Ploidy_primary_mean, dat$Ploidy_Met_mean, paired = TRUE)
cor.test(dat$Ploidy_primary_mean, dat$Ploidy_Met_mean, method = "spearman")
wilcox.test(dat$Purity_primary_mean, dat$Purity_Nodule_mean, paired = TRUE)
cor.test(dat$Purity_primary_mean, dat$Purity_Nodule_mean, method = "spearman")
wilcox.test(dat$Purity_primary_mean, dat$Purity_Lymph_mean, paired = TRUE)
cor.test(dat$Purity_primary_mean, dat$Purity_Lymph_mean, method = "spearman")
wilcox.test(dat$Purity_primary_mean, dat$Purity_Met_mean, paired = TRUE)
cor.test(dat$Purity_primary_mean, dat$Purity_Met_mean, method = "spearman")
library(ggplot2)
library(ggprism)
sda <- rbind(data.frame(id=dat$id, Ploidy=dat$Ploidy_primary_mean, Purity=dat$Purity_primary_mean, Class='P'),
  data.frame(id=dat$id, Ploidy=dat$Ploidy_Nodule_mean, Purity=dat$Purity_Nodule_mean, Class='N'),
  data.frame(id=dat$id, Ploidy=dat$Ploidy_Lymph_mean, Purity=dat$Purity_Lymph_mean, Class='LN')
  )
sda <- subset(sda, !is.na(Ploidy))
setheme <- theme(axis.text.y = element_text(size=16, color = "black"),
		  axis.text.x = element_text(size=16, color = "black"),
		  axis.title.y = element_text(size=16, color = "black", face = "bold"),
		  axis.title.x = element_text(size=16, color = "black", face = "bold"),
		  plot.title = element_text(size=16, color = "black",hjust = 0.5,margin = margin(t = 1, r = 0, b = 0.5, l = 0,"cm"), , face = "bold")) +
	theme(axis.line = element_line(colour = 'black', size = 0.5)) +
	theme(plot.margin = margin(t=0.2, r=2, b=0.2, l=1, "cm"))
sda$Type <- ifelse(sda$Class=='P','P\n(n=11)',ifelse(sda$Class=='LN','LN\n(n=3)',ifelse(sda$Class=='N','A\n(n=11)',NA)))
sda$Type <- factor(sda$Type,levels=c('P\n(n=11)','A\n(n=11)','LN\n(n=3)'),labels=c('P\n(n=11)','A\n(n=11)','LN\n(n=3)'))
p0 <- ggplot(aes(x=Type,y=Ploidy),data=sda) +
      geom_boxplot(fill='#3F76B4',color='#212121',width = 0.6,alpha=0.2,outlier.shape = NA) +
	  geom_jitter(width=0.1,color='#212121',size=3.5,alpha=0.3) +
	  theme_classic() +
	  xlab('') +
	  ylab('Ploidy') +
	  ylim(1.5,3.5) +
	  setheme
wilcox.test(dat$Ploidy_primary_mean, dat$Ploidy_Nodule_mean, paired = TRUE)
wilcox.test(dat$Ploidy_primary_mean, dat$Ploidy_Lymph_mean, paired = TRUE)
p_vals <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~label,
  "P\n(n=11)",   "A\n(n=11)",     0.76, 3.1, "p.exprs.ital",
  "P\n(n=11)",   "LN\n(n=3)",     0.50, 3.4, "p.exprs.ital"
)
p1 = p0 + add_pvalue(p_vals, label = "P = {p.adj}", tip.length = 0.02, label.size = 5)
p1
ggsave('Ploidy_Primary.VS.Nodule_paired.pdf',p1, width=5 ,height= 3.6)
p0 <- ggplot(aes(x=Type,y=Purity),data=sda) +
      geom_boxplot(fill='#3F76B4',color='#212121',width = 0.6,alpha=0.2,outlier.shape = NA) +
	  geom_jitter(width=0.1,color='#212121',size=3.5,alpha=0.3) +
	  theme_classic() +
	  xlab('') +
	  ylab('Purity') +
	  ylim(0.15,1) +
	  setheme
wilcox.test(dat$Purity_primary_mean, dat$Purity_Nodule_mean, paired = TRUE)
wilcox.test(dat$Purity_primary_mean, dat$Purity_Lymph_mean, paired = TRUE)
p_vals <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~label,
  "P\n(n=11)",   "A\n(n=11)",     0.37, 0.8, "p.exprs.ital",
  "P\n(n=11)",   "LN\n(n=3)",     1.00, 0.95, "p.exprs.ital"
)
p1 = p0 + add_pvalue(p_vals, label = "P = {p.adj}", tip.length = 0.02, label.size = 5)
p1
ggsave('Purity_Primary.VS.Met_paired.pdf',p1, width=5 ,height= 3.6)
