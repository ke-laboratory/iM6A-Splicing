rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(ggpubr)

#### Human
## Scatter plot
# Start
Data <- read.csv("iM6A_Human_SecondToLastExon.csv")
Data$Dvalue <- Data$Deletion - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
Data <- unique(Data)
table(Data$type)

ggplot(Data,aes(x = WT,y = Deletion), ) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

Gene <- Data[, c("chrom","name","strand")]
Gene <- unique(Gene)

# <200nt
Data1 <- subset(Data, SecondToLastExonLength<200)
table(Data1$type)

ggplot(Data1, aes(x = WT,y = Deletion), ) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(0,200) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
Gene <- Data1[, c("chrom","name","strand")]
Gene <- unique(Gene)

# 200-400nt
Data1 <- subset(Data, SecondToLastExonLength<=400&SecondToLastExonLength>=200)
table(Data1$type)

ggplot(Data1, aes(x = WT,y = Deletion), ) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(0,400) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
Gene <- Data1[, c("chrom","name","strand")]
Gene <- unique(Gene)

# >=400nt
Data1 <- subset(Data, SecondToLastExonLength>=400)
Data1 <- subset(Data1, Distance<=500)
table(Data1$type)

Gene <- Data1[, c("chrom","name","strand")]
Gene <- unique(Gene)

ggplot(Data1, aes(x = WT,y = Deletion), ) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

Data <- subset(Data, Distance<=500)
table(Data$type)

# End
Data <- read.csv("iM6A_Human_SecondToLastExon.csv")
Data$Dvalue <- Data$Deletion - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
Data <- unique(Data)
table(Data$type)
Data$Distance <- -(Data$SecondToLastExonLength - Data$Distance)

ggplot(Data,aes(x = WT,y = Deletion), ) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(-500,0) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")


# <200nt
Data1 <- subset(Data, SecondToLastExonLength<200)
table(Data1$type)

ggplot(Data1, aes(x = WT,y = Deletion), ) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(-200,0) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

# 200-400nt
Data1 <- subset(Data, SecondToLastExonLength<=400&SecondToLastExonLength>=200)
table(Data1$type)

ggplot(Data1, aes(x = WT,y = Deletion), ) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(-400,0) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

# >=400nt
Data1 <- subset(Data, SecondToLastExonLength>=400)
Data1 <- subset(Data1, Distance>= -500)
table(Data1$type)

ggplot(Data1, aes(x = WT,y = Deletion), ) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(-500,0) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

Data <- subset(Data, Distance>= -500)
table(Data$type)

#### Mouse
## Scatter plot
# Start
Data <- read.csv("iM6A_Mouse_SecondToLastExon.csv")
Data$Dvalue <- Data$Deletion - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
Data <- unique(Data)
table(Data$type)

ggplot(Data,aes(x = WT,y = Deletion), ) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

Gene <- Data[, c("chrom","name","strand")]
Gene <- unique(Gene)

# <200nt
Data1 <- subset(Data, SecondToLastExonLength<200)
table(Data1$type)

ggplot(Data1, aes(x = WT,y = Deletion), ) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(0,200) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

Gene <- Data1[, c("chrom","name","strand")]
Gene <- unique(Gene)

# 200-400nt
Data1 <- subset(Data, SecondToLastExonLength<=400&SecondToLastExonLength>=200)
table(Data1$type)

ggplot(Data1, aes(x = WT,y = Deletion), ) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(0,400) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

Gene <- Data1[, c("chrom","name","strand")]
Gene <- unique(Gene)

# >=400nt
Data1 <- subset(Data, SecondToLastExonLength>=400)
Data1 <- subset(Data1, Distance<=500)
table(Data1$type)

Gene <- Data1[, c("chrom","name","strand")]
Gene <- unique(Gene)

ggplot(Data1, aes(x = WT,y = Deletion), ) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

Data <- subset(Data, Distance<=500)
table(Data$type)



# End
Data <- read.csv("iM6A_Mouse_SecondToLastExon.csv")
Data$Dvalue <- Data$Deletion - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
Data <- unique(Data)
table(Data$type)
Data$Distance <- -(Data$SecondToLastExonLength - Data$Distance)

ggplot(Data,aes(x = WT,y = Deletion), ) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(-500,0) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")


# <200nt
Data1 <- subset(Data, SecondToLastExonLength<200)
table(Data1$type)

ggplot(Data1, aes(x = WT,y = Deletion), ) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(-200,0) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

# 200-400nt
Data1 <- subset(Data, SecondToLastExonLength<=400&SecondToLastExonLength>=200)
table(Data1$type)

ggplot(Data1, aes(x = WT,y = Deletion), ) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(-400,0) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

# >=400nt
Data1 <- subset(Data, SecondToLastExonLength>=400)
Data1 <- subset(Data1, Distance>= -500)
table(Data1$type)

ggplot(Data1, aes(x = WT,y = Deletion), ) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

ggplot()+geom_point(data=Data1, aes(x=Distance, y=Dvalue, group=Distance, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(-500,0) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

Data <- subset(Data, Distance>= -500)
table(Data$type)