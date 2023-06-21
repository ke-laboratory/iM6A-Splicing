rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(ggpubr)
library(ggExtra) 

#### Human
## Scatter plot
Cluster <- read.csv("Human_Cluster_LastExon.csv")
LastIntron <- read.csv("Human_LastIntron.csv")
Positive <- subset(LastIntron, strand=="+")
Negative <- subset(LastIntron, strand=="-")
Positive$LastIntronLength = Positive$LastExonStart - Positive$SecondToLastExonEnd
Negative$LastIntronLength = Negative$SecondToLastExonStart - Negative$LastExonEnd
LastIntron <- rbind(Positive, Negative)
LastIntron <- subset(LastIntron, LastIntronLength>5)
LastIntron <- merge(LastIntron, Cluster, by="name")

Data <- read.csv("iM6A_humanWhistleRAC_LastExon_LastIntronDeletion.csv")
Data <- merge(Data, LastIntron, by=c("chrom","name","strand"))

Gene <- Data[, c("chrom","name","strand")]
Gene <- unique(Gene)
Data <- merge(Data, Gene, by=c("chrom","name","strand"))

Data$Dvalue <- Data$Mutation-Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"

Positive <- subset(Data, strand=="+")
Negative <- subset(Data, strand=="-")
Positive$Distance <- Positive$Start-Positive$LastExonStart + 1
Negative$Distance <- Negative$LastExonEnd-Negative$End + 1
Data <- rbind(Positive, Negative)


# <=200nt
Data1 <- subset(Data, LastExonLength<200)

Gene <- Data1[, c("chrom","name","strand")]
Gene <- unique(Gene)
Data1 <- merge(Data1, Gene, by=c("chrom","name","strand"))


ggplot(Data1, aes(x = WT,y = Mutation)) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
  

ggplot()+geom_point(data=Data1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=1) + ylim(-1,1) + xlim(0,200) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
table(Data1$type)


# 200-400nt
Data1 <- subset(Data, LastExonLength>=200 & LastExonLength<=400)

Gene <- Data1[, c("chrom","name","strand")]
Gene <- unique(Gene)
set.seed(1)
Gene <- Gene %>% sample_n(1000, replace =FALSE)
Data1 <- merge(Data1, Gene, by=c("chrom","name","strand"))

ggplot(Data1, aes(x = WT,y = Mutation)) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none") 

ggplot()+geom_point(data=Data1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.5) + ylim(-1,1) + xlim(0,400) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
table(Data1$type)

# >400nt
Data1 <- subset(Data, LastExonLength>400)

Gene <- Data1[, c("chrom","name","strand")]
Gene <- unique(Gene)
set.seed(1)
Gene <- Gene %>% sample_n(1000, replace =FALSE)
Data1 <- merge(Data1, Gene, by=c("chrom","name","strand"))

ggplot(Data1, aes(x = WT,y = Mutation)) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.5) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
Data1 <- subset(Data1, Distance<=500)
table(Data1$type)


#### Mouse
## Scatter plot
Cluster <- read.csv("Mouse_Cluster_LastExon.csv")
LastIntron <- read.csv("Mouse_LastIntron.csv")
Positive <- subset(LastIntron, strand=="+")
Negative <- subset(LastIntron, strand=="-")
Positive$LastIntronLength = Positive$LastExonStart - Positive$SecondToLastExonEnd
Negative$LastIntronLength = Negative$SecondToLastExonStart - Negative$LastExonEnd
LastIntron <- rbind(Positive, Negative)
LastIntron <- subset(LastIntron, LastIntronLength>5)
LastIntron <- merge(LastIntron, Cluster, by="name")

Data <- read.csv("iM6A_LastExon_LastIntronDeletion.csv")
Data <- merge(Data, LastIntron, by=c("chrom","name","strand"))

Gene <- Data[, c("chrom","name","strand")]
Gene <- unique(Gene)
Data <- merge(Data, Gene, by=c("chrom","name","strand"))

Data$Dvalue <- Data$Mutation-Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"

Positive <- subset(Data, strand=="+")
Negative <- subset(Data, strand=="-")
Positive$Distance <- Positive$Start-Positive$LastExonStart + 1
Negative$Distance <- Negative$LastExonEnd-Negative$End + 1
Data <- rbind(Positive, Negative)


# <=200nt
Data1 <- subset(Data, LastExonLength<200)

Gene <- Data1[, c("chrom","name","strand")]
Gene <- unique(Gene)
Data1 <- merge(Data1, Gene, by=c("chrom","name","strand"))


ggplot(Data1, aes(x = WT,y = Mutation)) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")


ggplot()+geom_point(data=Data1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=1) + ylim(-1,1) + xlim(0,200) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
table(Data1$type)


# 200-400nt
Data1 <- subset(Data, LastExonLength>=200 & LastExonLength<=400)

Gene <- Data1[, c("chrom","name","strand")]
Gene <- unique(Gene)
set.seed(1)
Gene <- Gene %>% sample_n(1000, replace =FALSE)
Data1 <- merge(Data1, Gene, by=c("chrom","name","strand"))

ggplot(Data1, aes(x = WT,y = Mutation)) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none") 

ggplot()+geom_point(data=Data1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.5) + ylim(-1,1) + xlim(0,400) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
table(Data1$type)

# >400nt
Data1 <- subset(Data, LastExonLength>400)

Gene <- Data1[, c("chrom","name","strand")]
Gene <- unique(Gene)
set.seed(1)
Gene <- Gene %>% sample_n(1000, replace =FALSE)
Data1 <- merge(Data1, Gene, by=c("chrom","name","strand"))

ggplot(Data1, aes(x = WT,y = Mutation)) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.5) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
Data1 <- subset(Data1, Distance<=500)
table(Data1$type)
