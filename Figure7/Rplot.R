####################--------------------
rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(ggpubr)

## AA
Data <- read.csv("AA_InternalExon.csv")
Data <- Data[c(1,2,3,18,19,20),]
Data$AA <- factor(Data$AA, c("D", "N", "T", "A", "S", "P"))
ggplot(data=Data, aes(x=AA, y=OR)) + ylim(0,1.6) +
  geom_bar(stat="identity", color="blue", fill="white")

Data <- read.csv("AA_LastExon.csv")
Data <- Data[c(1,2,3,18,19,20),]
Data$AA <- factor(Data$AA, c("D", "N", "T", "S", "P", "A"))
ggplot(data=Data, aes(x=AA, y=OR)) + ylim(0,1.6) +
  geom_bar(stat="identity", color="blue", fill="white")


Data1 <- read.csv("AA_InternalExon.csv")
Data2 <- read.csv("AA_LastExon.csv")
Data1 <- Data1[, c(1,5)]
Data2 <- Data2[, c(1,5)]
colnames(Data1)[2] <- "IE"
colnames(Data2)[2] <- "LE"

Data <- merge(Data1, Data2, by="AA")
Data$IE <- log2(Data$IE)
Data$LE <- log2(Data$LE)
x = as.numeric(Data$IE)
y = as.numeric(Data$LE)
df1 = as.data.frame(cbind(x,y))
corT = cor.test(x,y,method="pearson")
cor = corT$estimate
pValue = corT$p.value
ggplot(Data, aes(x=IE, y=LE)) + geom_point(size=0.5) + xlim(-0.6,0.6)  + ylim(-0.6,0.6) + 
  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'pearson', aes(x=IE, y=LE))


####################--------------------
rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(ggpubr)

## Codon
Data <- read.csv("Codon_InternalExon.csv")
Data <- Data[c(1,2,3,4,5,57,58,59,60,61),]
Data$Codon <- factor(Data$Codon, c("GAC", "ACT", "GGA", "AAC", "ATA", "CGC", "GCG", "CCC", "TCG", "CCG"))
ggplot(data=Data, aes(x=Codon, y=OR)) + ylim(0,2) +
  geom_bar(stat="identity", color="blue", fill="white")

Data <- read.csv("Codon_LastExon.csv")
Data <- Data[c(1,2,3,4,5,57,58,59,60,61),]
Data$Codon <- factor(Data$Codon, c("GAC", "ACT", "AAC", "AGA", "GGA", "CGC", "GCC", "TCC", "TCG", "CCG"))
ggplot(data=Data, aes(x=Codon, y=OR)) + ylim(0,2) +
  geom_bar(stat="identity", color="blue", fill="white")


Data1 <- read.csv("Codon_InternalExon.csv")
Data2 <- read.csv("Codon_LastExon.csv")
Data1 <- Data1[, c(1,7)]
Data2 <- Data2[, c(1,7)]
colnames(Data1)[2] <- "IE"
colnames(Data2)[2] <- "LE"

Data <- merge(Data1, Data2, by="Codon")
Data$IE <- log2(Data$IE)
Data$LE <- log2(Data$LE)
x = as.numeric(Data$IE)
y = as.numeric(Data$LE)
df1 = as.data.frame(cbind(x,y))
corT = cor.test(x,y,method="pearson")
cor = corT$estimate
pValue = corT$p.value
ggplot(Data, aes(x=IE, y=LE)) + geom_point(size=0.5) + xlim(-1,1)  + ylim(-1,1) + 
  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'pearson', aes(x=IE, y=LE))

####################--------------------
rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(ggpubr)
Data <- read.csv("Codon_InternalExon.csv")
Data <- subset(Data, Codon=="GAC" | Codon=="GAT")
fisher.test( Data[,c(2,3)])
Data$Codon <- factor(Data$Codon, c("GAC","GAT"))
ggplot(data=Data, aes(x=Codon, y=OR)) + ylim(0,2) +
  geom_bar(stat="identity", color="blue", fill="white")

Data <- read.csv("Codon_InternalExon.csv")
Data <- subset(Data, Codon=="AAC" | Codon=="AAT")
fisher.test( Data[,c(2,3)])
Data$Codon <- factor(Data$Codon, c("AAC","AAT"))
ggplot(data=Data, aes(x=Codon, y=OR)) + ylim(0,2) +
  geom_bar(stat="identity", color="blue", fill="white")


###########-------------
rm(list=ls())
options(stringsAsFactors = F)
library(ggpubr)
library(plyr)

# C1 VS C2
Data <- read.csv("Cluster_Length.csv")
Data <- subset(Data, ExonLength>=50)
C1 <- subset(Data, Cluster=="C1")
C2 <- subset(Data, Cluster=="C2")

ggplot(Data, aes(x=ExonLength , color=Cluster))  + 
  geom_density(alpha=0.6) + scale_fill_grey() + theme_classic() + xlim(0,500)
sample1 <- C1$ExonLength
sample2 <- C2$ExonLength
ks.test(sample1, sample2)


ggplot(Data, aes(x=IntronLength1 , color=Cluster))  + 
  geom_density(alpha=0.6) + scale_fill_grey() + theme_classic() + xlim(0,5000)
sample1 <- C1$IntronLength1
sample2 <- C2$IntronLength1
ks.test(sample1, sample2)


ggplot(Data, aes(x=IntronLength2 , color=Cluster)) + xlim(0,5000) + 
  geom_density(alpha=0.6) + scale_fill_grey() + theme_classic()
sample1 <- C1$IntronLength2
sample2 <- C2$IntronLength2
ks.test(sample1, sample2)


