rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(ggpubr)
library(ggExtra) 

#### Human
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
Data <- unique(Data)

C1 <- subset(Data, Cluster=="C1")
C2 <- subset(Data, Cluster=="C2")
C1 <- C1[, c(2,26)]
C2 <- C2[, c(2,26)]
C1 <- unique(C1)
C2 <- unique(C2)
set.seed(1)
C1 <- C1 %>% sample_n(2000, replace =FALSE)
set.seed(1)
C2 <- C2 %>% sample_n(2000, replace =FALSE)

Cluster <- rbind(C1, C2)
Data <- merge(Data, Cluster, by=c("name","Cluster"))

Data$Dvalue <- Data$Mutation-Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"

Positive <- subset(Data, strand=="+")
Negative <- subset(Data, strand=="-")
Positive$Distance <- Positive$Start-Positive$LastExonStart + 1
Negative$Distance <- Negative$LastExonEnd-Negative$End + 1
Data <- rbind(Positive, Negative)

# C1
C1 <- subset(Data, Cluster=="C1")
table(C1$type)
ggplot(C1,aes(x = WT,y = Mutation)) + xlim(0,1)  + ylim(0,1) + 
  geom_point(aes(color = type), size=1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) +
  geom_abline(slope = 1,intercept = 0.1,lty = 'dashed',size = 0.8) +
  geom_abline(slope = 1,intercept = -0.1,lty = 'dashed',size = 0.8) + theme(legend.position = "none")

ggplot()+geom_point(data=C1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

# C2
C2 <- subset(Data, Cluster=="C2")
table(C2$type)
ggplot(C2,aes(x = WT,y = Mutation)) + xlim(0,1)  + ylim(0,1) + 
  geom_point(aes(color = type), size=1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) +
  geom_abline(slope = 1,intercept = 0.1,lty = 'dashed',size = 0.8) +
  geom_abline(slope = 1,intercept = -0.1,lty = 'dashed',size = 0.8) + theme(legend.position = "none")

ggplot()+geom_point(data=C2, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

#### Mouse
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
Data <- unique(Data)

C1 <- subset(Data, Cluster=="C1")
C2 <- subset(Data, Cluster=="C2")
C1 <- C1[, c(2,26)]
C2 <- C2[, c(2,26)]
C1 <- unique(C1)
C2 <- unique(C2)
set.seed(1)
C1 <- C1 %>% sample_n(2000, replace =FALSE)
set.seed(1)
C2 <- C2 %>% sample_n(2000, replace =FALSE)

Cluster <- rbind(C1, C2)
Data <- merge(Data, Cluster, by=c("name","Cluster"))

Data$Dvalue <- Data$Mutation-Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"

Positive <- subset(Data, strand=="+")
Negative <- subset(Data, strand=="-")
Positive$Distance <- Positive$Start-Positive$LastExonStart + 1
Negative$Distance <- Negative$LastExonEnd-Negative$End + 1
Data <- rbind(Positive, Negative)

# C1
C1 <- subset(Data, Cluster=="C1")
ggplot(C1,aes(x = WT,y = Mutation)) + xlim(0,1)  + ylim(0,1) + 
  geom_point(aes(color = type), size=1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) +
  geom_abline(slope = 1,intercept = 0.1,lty = 'dashed',size = 0.8) +
  geom_abline(slope = 1,intercept = -0.1,lty = 'dashed',size = 0.8) + theme(legend.position = "none")

ggplot()+geom_point(data=C1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

# C2
C2 <- subset(Data, Cluster=="C2")
table(C2$type)
ggplot(C2,aes(x = WT,y = Mutation)) + xlim(0,1)  + ylim(0,1) + 
  geom_point(aes(color = type), size=1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) +
  geom_abline(slope = 1,intercept = 0.1,lty = 'dashed',size = 0.8) +
  geom_abline(slope = 1,intercept = -0.1,lty = 'dashed',size = 0.8) + theme(legend.position = "none")

ggplot()+geom_point(data=C2, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

