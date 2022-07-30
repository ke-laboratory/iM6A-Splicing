###########-------------
rm(list=ls())
options(stringsAsFactors = F)
library(ggpubr)
library(dplyr)
library(plyr)

Data <- read.csv("Data.csv")
colnames(Data)[9] <- "FC"
Data$Bin[Data$exonCount==1] <-"Bin0"
Data$Bin[Data$exonCount==2] <-"Bin1"
Data$Bin[Data$exonCount==3] <-"Bin2"
Data$Bin[Data$exonCount==4] <-"Bin3"
Data$Bin[Data$exonCount==5] <-"Bin4"
Data$Bin[Data$exonCount==6] <-"Bin5"
Data$Bin[Data$exonCount==7] <-"Bin6"
Data$Bin[Data$exonCount==8] <-"Bin7"
Data$Bin[Data$exonCount==9] <-"Bin8"
Data$Bin[Data$exonCount>=10] <-"Bin9"
Data$M6ADensity <- Data$Count*100/Data$cDNALength
Data$ExonDensity <- Data$exonCount*100/Data$cDNALength
RAC <- read.csv("Count_RAC.csv")
Data <- merge(Data, RAC, by="name")
Data$Percent <- Data$Count/Data$RAC

## Scatter plot
x = as.numeric(Data$exonCount)
y = as.numeric(Data$Percent)
df1 = as.data.frame(cbind(x,y))
corT = cor.test(x,y,method="pearson")
cor = corT$estimate
pValue = corT$p.value
ggplot(Data, aes(x=exonCount, y=Percent)) + geom_point(size=0.5) + xlim(0,60)  + ylim(0,0.25) + 
  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'pearson', aes(x=exonCount, y=Percent))

x = as.numeric(Data$WT)
y = as.numeric(Data$Percent)
df1 = as.data.frame(cbind(x,y))
corT = cor.test(x,y,method="pearson")
cor = corT$estimate
pValue = corT$p.value
ggplot(Data, aes(x=WT, y=Percent)) + geom_point(size=0.5) + xlim(0,20)  + ylim(0,0.25) + 
  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'pearson', aes(x=WT, y=Percent))

x = as.numeric(Data$KO)
y = as.numeric(Data$Percent)
df1 = as.data.frame(cbind(x,y))
corT = cor.test(x,y,method="pearson")
cor = corT$estimate
pValue = corT$p.value
ggplot(Data, aes(x=KO, y=Percent)) + geom_point(size=0.5) + xlim(0,20)  + ylim(0,0.25) + 
  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'pearson', aes(x=KO, y=Percent))

x = as.numeric(Data$WT)
y = as.numeric(Data$exonCount)
df1 = as.data.frame(cbind(x,y))
corT = cor.test(x,y,method="pearson")
cor = corT$estimate
pValue = corT$p.value
ggplot(Data, aes(x=WT, y=exonCount)) + geom_point(size=0.5) + xlim(0,20)  + ylim(0,60) + 
  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'pearson', aes(x=WT, y=exonCount))

x = as.numeric(Data$KO)
y = as.numeric(Data$exonCount)
df1 = as.data.frame(cbind(x,y))
corT = cor.test(x,y,method="pearson")
cor = corT$estimate
pValue = corT$p.value
ggplot(Data, aes(x=KO, y=exonCount)) + geom_point(size=0.5) + xlim(0,20)  + ylim(0,60) + 
  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'pearson', aes(x=KO, y=exonCount))

## Box plot
ggplot(Data, aes(x=Bin, y=WT, fill=Bin)) + geom_boxplot(outlier.colour=NA) + geom_jitter(shape=16, position=position_jitter(0.2), size=1) + ylim(0,15) + theme(legend.position = "none")
ggplot(Data, aes(x=Bin, y=KO, fill=Bin)) + geom_boxplot(outlier.colour=NA) + geom_jitter(shape=16, position=position_jitter(0.2), size=1) + ylim(0,15) + theme(legend.position = "none")
ggplot(Data, aes(x=Bin, y=M6ADensity, fill=Bin)) + geom_boxplot(outlier.colour=NA) + geom_jitter(shape=16, position=position_jitter(0.2), size=1) + ylim(0,1) + theme(legend.position = "none")
ggplot(Data, aes(x=Bin, y=Percent, fill=Bin)) + geom_boxplot(outlier.colour=NA) + geom_jitter(shape=16, position=position_jitter(0.2), size=1) + ylim(0,0.25) + theme(legend.position = "none")




###########-------------
rm(list=ls())
options(stringsAsFactors = F)
library(ggpubr)
library(dplyr)
library(plyr)

Data <- read.csv("Data_Single_VS_Multiple.csv")
colnames(Data)[10] <- "FC"
Data$Mark[Data$exonCount<=1] <- "A"
Data$Mark[Data$exonCount>1] <- "B"
A <- subset(Data, Mark=="A")
B <- subset(Data, Mark=="B")


# Half-life of WT
x=colnames(Data)[11]
y=colnames(Data)[9]
group=levels(factor(Data$Mark))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Data, x="Mark", y="WT", fill = "Mark", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons)
ggplot(Data, aes(x=Mark, y=WT, fill=Mark)) + geom_boxplot(outlier.colour=NA)+ ylim(0,15) + theme(legend.position = "none")

sample1 <- A$WT
sample2 <- B$WT
ks.test(sample1, sample2)
group <- c(rep("sample1", length(sample1)), rep("sample2", length(sample2)))
dat <- data.frame(KSD = c(sample1,sample2), group = group)
cdf1 <- ecdf(sample1) 
cdf2 <- ecdf(sample2) 
ggplot(dat, aes(x = KSD, group = group, color = group)) + stat_ecdf(aes(colour=group)) + xlim(0,40) + theme(legend.position = "none")

# Half-life of KO
x=colnames(Data)[10]
y=colnames(Data)[8]
group=levels(factor(Data$Mark))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Data, x="Mark", y="KO", fill = "Mark", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons)
ggplot(Data, aes(x=Mark, y=KO, fill=Mark)) + geom_boxplot(outlier.colour=NA)+ ylim(0,20) + theme(legend.position = "none")

sample1 <- A$KO
sample2 <- B$KO
ks.test(sample1, sample2)
group <- c(rep("sample1", length(sample1)), rep("sample2", length(sample2)))
dat <- data.frame(KSD = c(sample1,sample2), group = group)
cdf1 <- ecdf(sample1) 
cdf2 <- ecdf(sample2) 
ggplot(dat, aes(x = KSD, group = group, color = group)) + stat_ecdf(aes(colour=group)) + xlim(0,40) + theme(legend.position = "none")

# FC
x=colnames(Data)[11]
y=colnames(Data)[10]
group=levels(factor(Data$Mark))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Data, x="Mark", y="FC", fill = "Mark", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons)
ggplot(Data, aes(x=Mark, y=FC, fill=Mark)) + geom_boxplot(outlier.colour=NA) + theme(legend.position = "none")

sample1 <- A$FC
sample2 <- B$FC
ks.test(sample1, sample2)
group <- c(rep("sample1", length(sample1)), rep("sample2", length(sample2)))
dat <- data.frame(KSD = c(sample1,sample2), group = group)
cdf1 <- ecdf(sample1) 
cdf2 <- ecdf(sample2) 
ggplot(dat, aes(x = KSD, group = group, color = group)) + stat_ecdf(aes(colour=group)) + xlim(-3,3) + theme(legend.position = "none")

# m6A
x=colnames(Data)[11]
y=colnames(Data)[6]
group=levels(factor(Data$Mark))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Data, x="Mark", y="Count", fill = "Mark", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons)
ggplot(Data, aes(x=Mark, y=Count, fill=Mark)) + geom_boxplot(outlier.colour=NA) + ylim(0,15) + theme(legend.position = "none")

sample1 <- A$Count
sample2 <- B$Count
ks.test(sample1, sample2)
group <- c(rep("sample1", length(sample1)), rep("sample2", length(sample2)))
dat <- data.frame(KSD = c(sample1,sample2), group = group)
cdf1 <- ecdf(sample1) 
cdf2 <- ecdf(sample2) 
ggplot(dat, aes(x = KSD, group = group, color = group)) + stat_ecdf(aes(colour=group)) + theme(legend.position = "none")

### Less VS More
rm(list=ls())
options(stringsAsFactors = F)
library(ggpubr)
library(dplyr)
library(plyr)

Data <- read.csv("Data_Less_VS_More.csv")
colnames(Data)[10] <- "FC"
Data$Mark[Data$exonCount<=6] <- "A"
Data$Mark[Data$exonCount>6] <- "B"
A <- subset(Data, Mark=="A")
B <- subset(Data, Mark=="B")


# Half-life for WT
x=colnames(Data)[11]
y=colnames(Data)[9]
group=levels(factor(Data$Mark))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Data, x="Mark", y="WT", fill = "Mark", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons)
ggplot(Data, aes(x=Mark, y=WT, fill=Mark)) + geom_boxplot(outlier.colour=NA)+ ylim(0,15) + theme(legend.position = "none")

sample1 <- A$WT
sample2 <- B$WT
ks.test(sample1, sample2)
group <- c(rep("sample1", length(sample1)), rep("sample2", length(sample2)))
dat <- data.frame(KSD = c(sample1,sample2), group = group)
cdf1 <- ecdf(sample1) 
cdf2 <- ecdf(sample2) 
ggplot(dat, aes(x = KSD, group = group, color = group)) + stat_ecdf(aes(colour=group)) + xlim(0,40)+ theme(legend.position = "none")

# Half-life for KO
x=colnames(Data)[11]
y=colnames(Data)[8]
group=levels(factor(Data$Mark))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Data, x="Mark", y="KO", fill = "Mark", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons)
ggplot(Data, aes(x=Mark, y=KO, fill=Mark)) + geom_boxplot(outlier.colour=NA)+ ylim(0,20)+ theme(legend.position = "none")

sample1 <- A$KO
sample2 <- B$KO
ks.test(sample1, sample2)
group <- c(rep("sample1", length(sample1)), rep("sample2", length(sample2)))
dat <- data.frame(KSD = c(sample1,sample2), group = group)
cdf1 <- ecdf(sample1) 
cdf2 <- ecdf(sample2) 
ggplot(dat, aes(x = KSD, group = group, color = group)) + stat_ecdf(aes(colour=group)) + xlim(0,40)+ theme(legend.position = "none")

# FC
x=colnames(Data)[11]
y=colnames(Data)[10]
group=levels(factor(Data$Mark))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Data, x="Mark", y="FC", fill = "Mark", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons)
ggplot(Data, aes(x=Mark, y=FC, fill=Mark)) + geom_boxplot(outlier.colour=NA) + ylim(-2,4)+ theme(legend.position = "none")

sample1 <- A$FC
sample2 <- B$FC
ks.test(sample1, sample2)
group <- c(rep("sample1", length(sample1)), rep("sample2", length(sample2)))
dat <- data.frame(KSD = c(sample1,sample2), group = group)
cdf1 <- ecdf(sample1) 
cdf2 <- ecdf(sample2) 
ggplot(dat, aes(x = KSD, group = group, color = group)) + stat_ecdf(aes(colour=group)) + xlim(-2,4)+ theme(legend.position = "none")

# m6A
x=colnames(Data)[11]
y=colnames(Data)[6]
group=levels(factor(Data$Mark))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Data, x="Mark", y="Count", fill = "Mark", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons)
ggplot(Data, aes(x=Mark, y=Count, fill=Mark)) + geom_boxplot(outlier.colour=NA) + ylim(0,10)+ theme(legend.position = "none")

sample1 <- A$Count
sample2 <- B$Count
ks.test(sample1, sample2)
group <- c(rep("sample1", length(sample1)), rep("sample2", length(sample2)))
dat <- data.frame(KSD = c(sample1,sample2), group = group)
cdf1 <- ecdf(sample1) 
cdf2 <- ecdf(sample2) 
ggplot(dat, aes(x = KSD, group = group, color = group)) + stat_ecdf(aes(colour=group)) + xlim(0,20) + theme(legend.position = "none")


### C1 VS C2
rm(list=ls())
options(stringsAsFactors = F)
library(ggpubr)
library(dplyr)
library(plyr)

Data <- read.csv("Data_C1_VS_C2.csv")
colnames(Data)[16] <- "FC"

Data$Mark[Data$C1>0] <- "A"
Data$Mark[Data$C1==0] <- "B"
Data$M6ADensity <- Data$Count*100/Data$cDNALength
Data$ExonDensity <- Data$exonCount*100/Data$cDNALength
Data$Percent <- Data$Count/Data$RAC

A <- subset(Data, Mark=="A")
B <- subset(Data, Mark=="B")

# Half-life for m6A/RAC percent
x=colnames(Data)[18]
y=colnames(Data)[21]
group=levels(factor(Data$Mark))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Data, x="Mark", y="Percent", fill = "Mark", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons)
ggplot(Data, aes(x=Mark, y=Percent, fill=Mark)) + geom_boxplot(outlier.colour=NA)+ ylim(0,0.1) + theme(legend.position = "none")

sample1 <- A$Percent
sample2 <- B$Percent
ks.test(sample1, sample2)
group <- c(rep("sample1", length(sample1)), rep("sample2", length(sample2)))
dat <- data.frame(KSD = c(sample1,sample2), group = group)
cdf1 <- ecdf(sample1) 
cdf2 <- ecdf(sample2) 
ggplot(dat, aes(x = KSD, group = group, color = group)) + stat_ecdf(aes(colour=group)) + xlim(0,0.2) + theme(legend.position = "none")


# Half-life for WT
x=colnames(Data)[17]
y=colnames(Data)[15]
group=levels(factor(Data$Mark))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Data, x="Mark", y="WT", fill = "Mark", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons)
ggplot(Data, aes(x=Mark, y=WT, fill=Mark)) + geom_boxplot(outlier.colour=NA)+ ylim(0,15) + theme(legend.position = "none")

sample1 <- A$WT
sample2 <- B$WT
ks.test(sample1, sample2)
group <- c(rep("sample1", length(sample1)), rep("sample2", length(sample2)))
dat <- data.frame(KSD = c(sample1,sample2), group = group)
cdf1 <- ecdf(sample1) 
cdf2 <- ecdf(sample2) 
ggplot(dat, aes(x = KSD, group = group, color = group)) + stat_ecdf(aes(colour=group)) + xlim(0,20) + theme(legend.position = "none")

# Half-life for KO
x=colnames(Data)[17]
y=colnames(Data)[14]
group=levels(factor(Data$Mark))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Data, x="Mark", y="KO", fill = "Mark", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons)
ggplot(Data, aes(x=Mark, y=KO, fill=Mark)) + geom_boxplot(outlier.colour=NA)+ ylim(0,15) + theme(legend.position = "none")

sample1 <- A$KO
sample2 <- B$KO
ks.test(sample1, sample2)
group <- c(rep("sample1", length(sample1)), rep("sample2", length(sample2)))
dat <- data.frame(KSD = c(sample1,sample2), group = group)
cdf1 <- ecdf(sample1) 
cdf2 <- ecdf(sample2) 
ggplot(dat, aes(x = KSD, group = group, color = group)) + stat_ecdf(aes(colour=group)) + xlim(0,20) + theme(legend.position = "none")

# FC
x=colnames(Data)[17]
y=colnames(Data)[16]
group=levels(factor(Data$Mark))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Data, x="Mark", y="FC", fill = "Mark", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons)
ggplot(Data, aes(x=Mark, y=FC, fill=Mark)) + geom_boxplot(outlier.colour=NA) + ylim(-2,3) + theme(legend.position = "none")

sample1 <- A$FC
sample2 <- B$FC
ks.test(sample1, sample2)
group <- c(rep("sample1", length(sample1)), rep("sample2", length(sample2)))
dat <- data.frame(KSD = c(sample1,sample2), group = group)
cdf1 <- ecdf(sample1) 
cdf2 <- ecdf(sample2) 
ggplot(dat, aes(x = KSD, group = group, color = group)) + stat_ecdf(aes(colour=group)) + xlim(-3,3) + theme(legend.position = "none")

# m6A
x=colnames(Data)[17]
y=colnames(Data)[13]
group=levels(factor(Data$Mark))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Data, x="Mark", y="Count", fill = "Mark", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons)
ggplot(Data, aes(x=Mark, y=Count, fill=Mark)) + geom_boxplot(outlier.colour=NA) + ylim(0,10) + theme(legend.position = "none")

sample1 <- A$Count
sample2 <- B$Count
ks.test(sample1, sample2)
group <- c(rep("sample1", length(sample1)), rep("sample2", length(sample2)))
dat <- data.frame(KSD = c(sample1,sample2), group = group)
cdf1 <- ecdf(sample1) 
cdf2 <- ecdf(sample2) 
ggplot(dat, aes(x = KSD, group = group, color = group)) + stat_ecdf(aes(colour=group)) + theme(legend.position = "none")













