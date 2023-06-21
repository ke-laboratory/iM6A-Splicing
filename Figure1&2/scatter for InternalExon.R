rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(ggpubr)

#### Human
## Scatter plot
# Start
Data <- read.csv("iM6A_Human_IL_Start.csv")
Data <- Data[, c("chrom","name","strand","m6AStart","m6AEnd","WT","IL","m6APosition","ExonLength")]
Data <- unique(Data)

Gene <- Data[, c("chrom","name","strand")]
Gene <- unique(Gene)
set.seed(1)
Gene <- Gene %>% sample_n(1000, replace =FALSE)

Data <- merge(Data, Gene, by=c("chrom","name","strand"))

Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
table(Data$type)

ggplot(Data,aes(x = WT,y = IL)) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none") +
  geom_abline(slope = 1,intercept = 0.1,lty = 'dashed',size = 0.8) +
  geom_abline(slope = 1,intercept = -0.1,lty = 'dashed',size = 0.8)

ggplot()+geom_point(data=Data, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
Data <- subset(Data, m6APosition<=500)
table(Data$type)

# End
Data$m6APosition <- -(Data$ExonLength - Data$m6APosition)
ggplot()+geom_point(data=Data, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(-500,0) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
Data <- subset(Data, m6APosition>=-500)
table(Data$type)

## <200nt 
# Start
Data <- read.csv("iM6A_Human_IL_Start.csv")
Data <- Data[, c("chrom","name","strand","m6AStart","m6AEnd","WT","IL","m6APosition","ExonLength")]
Data <- unique(Data)
Data <- subset(Data, ExonLength<200)

Gene <- Data[, c("chrom","name","strand")]
Gene <- unique(Gene)
set.seed(1)
Gene <- Gene %>% sample_n(1000, replace =FALSE)

Data <- merge(Data, Gene, by=c("chrom","name","strand"))

Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
table(Data$type)

ggplot(Data,aes(x = WT,y = IL)) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none") +
  geom_abline(slope = 1,intercept = 0.1,lty = 'dashed',size = 0.8) +
  geom_abline(slope = 1,intercept = -0.1,lty = 'dashed',size = 0.8)

ggplot()+geom_point(data=Data, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(0,200) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

# End
Data$m6APosition <- -(Data$ExonLength - Data$m6APosition)
ggplot()+geom_point(data=Data, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(-200,0) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

## 200-400nt 
# Start
Data <- read.csv("iM6A_Human_IL_Start.csv")
Data <- Data[, c("chrom","name","strand","m6AStart","m6AEnd","WT","IL","m6APosition","ExonLength")]
Data <- unique(Data)
Data <- subset(Data, ExonLength>=200&ExonLength<=400)

Gene <- Data[, c("chrom","name","strand")]
Gene <- unique(Gene)
set.seed(1)
Gene <- Gene %>% sample_n(1000, replace =FALSE)

Data <- merge(Data, Gene, by=c("chrom","name","strand"))

Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
table(Data$type)

ggplot(Data,aes(x = WT,y = IL)) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none") +
  geom_abline(slope = 1,intercept = 0.1,lty = 'dashed',size = 0.8) +
  geom_abline(slope = 1,intercept = -0.1,lty = 'dashed',size = 0.8)

ggplot()+geom_point(data=Data, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(0,400) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

# End
Data$m6APosition <- -(Data$ExonLength - Data$m6APosition)
ggplot()+geom_point(data=Data, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(-400,0) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

## >=400nt 
# Start
Data <- read.csv("iM6A_Human_IL_Start.csv")
Data <- Data[, c("chrom","name","strand","m6AStart","m6AEnd","WT","IL","m6APosition","ExonLength")]
Data <- unique(Data)
Data <- subset(Data, ExonLength>=400)

Gene <- Data[, c("chrom","name","strand")]
Gene <- unique(Gene)
set.seed(1)
Gene <- Gene %>% sample_n(2000, replace =FALSE)

Data <- merge(Data, Gene, by=c("chrom","name","strand"))

Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
table(Data$type)
  
ggplot(Data,aes(x = WT,y = IL)) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none") +
  geom_abline(slope = 1,intercept = 0.1,lty = 'dashed',size = 0.8) +
  geom_abline(slope = 1,intercept = -0.1,lty = 'dashed',size = 0.8)

ggplot()+geom_point(data=Data, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
Data <- subset(Data, m6APosition<=500)
table(Data$type)

# End
Data$m6APosition <- -(Data$ExonLength - Data$m6APosition)
ggplot()+geom_point(data=Data, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(-500,0) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
Data <- subset(Data, m6APosition>= -500)
table(Data$type)

#### Mouse
## Scatter plot
# Start
Data <- read.csv("iM6A_Mouse_IL_Start.csv")
Data <- Data[, c("chrom","name","strand","m6AStart","m6AEnd","WT","IL","m6APosition","ExonLength")]
Data <- unique(Data)

Gene <- Data[, c("chrom","name","strand")]
Gene <- unique(Gene)
set.seed(1)
Gene <- Gene %>% sample_n(1000, replace =FALSE)

Data <- merge(Data, Gene, by=c("chrom","name","strand"))

Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
table(Data$type)

ggplot(Data,aes(x = WT,y = IL)) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none") +
  geom_abline(slope = 1,intercept = 0.1,lty = 'dashed',size = 0.8) +
  geom_abline(slope = 1,intercept = -0.1,lty = 'dashed',size = 0.8)

ggplot()+geom_point(data=Data, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
Data <- subset(Data, m6APosition<=500)
table(Data$type)

# End
Data$m6APosition <- -(Data$ExonLength - Data$m6APosition)
ggplot()+geom_point(data=Data, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(-500,0) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
Data <- subset(Data, m6APosition>=-500)
table(Data$type)


## <200nt 
# Start
Data <- read.csv("iM6A_Mouse_IL_Start.csv")
Data <- Data[, c("chrom","name","strand","m6AStart","m6AEnd","WT","IL","m6APosition","ExonLength")]
Data <- unique(Data)
Data <- subset(Data, ExonLength<200)

Gene <- Data[, c("chrom","name","strand")]
Gene <- unique(Gene)
set.seed(1)
Gene <- Gene %>% sample_n(1000, replace =FALSE)

Data <- merge(Data, Gene, by=c("chrom","name","strand"))

Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
table(Data$type)

ggplot(Data,aes(x = WT,y = IL)) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none") +
  geom_abline(slope = 1,intercept = 0.1,lty = 'dashed',size = 0.8) +
  geom_abline(slope = 1,intercept = -0.1,lty = 'dashed',size = 0.8)

ggplot()+geom_point(data=Data, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(0,200) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

# End
Data$m6APosition <- -(Data$ExonLength - Data$m6APosition)
ggplot()+geom_point(data=Data, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(-200,0) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

## 200-400nt 
# Start
Data <- read.csv("iM6A_Mouse_IL_Start.csv")
Data <- Data[, c("chrom","name","strand","m6AStart","m6AEnd","WT","IL","m6APosition","ExonLength")]
Data <- unique(Data)
Data <- subset(Data, ExonLength>=200&ExonLength<=400)

Gene <- Data[, c("chrom","name","strand")]
Gene <- unique(Gene)
set.seed(1)
Gene <- Gene %>% sample_n(1000, replace =FALSE)

Data <- merge(Data, Gene, by=c("chrom","name","strand"))

Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
table(Data$type)

ggplot(Data,aes(x = WT,y = IL)) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none") +
  geom_abline(slope = 1,intercept = 0.1,lty = 'dashed',size = 0.8) +
  geom_abline(slope = 1,intercept = -0.1,lty = 'dashed',size = 0.8)

ggplot()+geom_point(data=Data, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(0,400) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

# End
Data$m6APosition <- -(Data$ExonLength - Data$m6APosition)
ggplot()+geom_point(data=Data, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(-400,0) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")

## >=400nt 
# Start
Data <- read.csv("iM6A_Mouse_IL_Start.csv")
Data <- Data[, c("chrom","name","strand","m6AStart","m6AEnd","WT","IL","m6APosition","ExonLength")]
Data <- unique(Data)
Data <- subset(Data, ExonLength>=400)

Gene <- Data[, c("chrom","name","strand")]
Gene <- unique(Gene)
set.seed(1)
Gene <- Gene %>% sample_n(2000, replace =FALSE)

Data <- merge(Data, Gene, by=c("chrom","name","strand"))

Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
table(Data$type)

ggplot(Data,aes(x = WT,y = IL)) + xlim(0,1) + ylim(0,1) + 
  geom_point(aes(color = type), size=0.1) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none") +
  geom_abline(slope = 1,intercept = 0.1,lty = 'dashed',size = 0.8) +
  geom_abline(slope = 1,intercept = -0.1,lty = 'dashed',size = 0.8)

ggplot()+geom_point(data=Data, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
Data <- subset(Data, m6APosition<=500)
table(Data$type)

# End
Data$m6APosition <- -(Data$ExonLength - Data$m6APosition)
ggplot()+geom_point(data=Data, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.25) + ylim(-1,1) + xlim(-500,0) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
Data <- subset(Data, m6APosition>=-500)
table(Data$type)
