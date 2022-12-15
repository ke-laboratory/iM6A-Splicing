rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(ggpubr)
library(ggExtra) 

## Plot for last exon
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
set.seed(1)
Gene <- Gene %>% sample_n(4000, replace =FALSE)

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

ggplot(Data,aes(x = WT,y = Mutation)) + xlim(0,1)  + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#3E7C17')) +
  geom_abline(slope = 1,intercept = 0.1,lty = 'dashed',size = 0.8) +
  geom_abline(slope = 1,intercept = -0.1,lty = 'dashed',size = 0.8) + theme(legend.position = "none")

ggplot()+geom_point(data=Data, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(0,1000) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#3E7C17')) + theme(legend.position = "none")


## Plot for second-to-last exon
# Start
Data <- read.csv("iM6A_Mouse_SecondToLastExon.csv")
Data$Dvalue <- Data$Deletion - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
Data <- unique(Data)

ggplot(Data,aes(x = WT,y = Deletion), ) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#3E7C17')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(0,1000) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#3E7C17')) + theme(legend.position = "none")

# End
Data <- read.csv("iM6A_Mouse_SecondToLastExon.csv")
Data$Dvalue <- Data$Deletion - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
Data <- unique(Data)
Data$Distance <- Data$SecondToLastExonLength - Data$Distance

ggplot(Data,aes(x = WT,y = Deletion), ) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#3E7C17')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(0,1000) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#3E7C17')) + theme(legend.position = "none")

## Plot for internal exon
# Start
Data <- read.csv("iM6A_Mouse_IL_Start.csv")
Data <- Data[, c("chrom","name","strand","m6AStart","m6AEnd","WT","IL","m6APosition","ExonLength")]
Data <- unique(Data)
Gene <- Data[, c("chrom","name","strand")]
Gene <- unique(Gene)
set.seed(1)
Gene <- Gene %>% sample_n(1500, replace =FALSE)
Data <- merge(Data, Gene, by=c("chrom","name","strand"))
Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"

ggplot(Data,aes(x = WT,y = IL)) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#3E7C17')) + theme(legend.position = "none") +
  geom_abline(slope = 1,intercept = 0.1,lty = 'dashed',size = 0.8) +
  geom_abline(slope = 1,intercept = -0.1,lty = 'dashed',size = 0.8)

ggplot()+geom_point(data=Data, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(0,1000) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#3E7C17')) + theme(legend.position = "none")

# End
Data$m6APosition <- Data$ExonLength - Data$m6APosition
ggplot()+geom_point(data=Data, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(0,1000) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#3E7C17')) + theme(legend.position = "none")








