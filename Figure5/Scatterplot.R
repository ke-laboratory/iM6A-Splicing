rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(ggpubr)

#### Mouse
## Scatter plot
# Start
Data <- read.csv("iM6A_Mouse_IL_Start.csv")
Data <- Data[, c("chrom","name","strand","ID","m6AStart","m6AEnd","WT","IL","m6APosition")]
Data <- unique(Data)

Cluster <- read.csv("Mouse_Cluster_InternalExon.csv")
C1 <- subset(Cluster, Cluster=="C1")
C2 <- subset(Cluster, Cluster=="C2")
set.seed(1)
C1 <- C1 %>% sample_n(3000, replace =FALSE)
set.seed(1)
C2 <- C2 %>% sample_n(3000, replace =FALSE)
Cluster <- rbind(C1,C2)

Data <- merge(Data, Cluster, by="ID")
Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"

# C1
C1 <- subset(Data, Cluster=="C1")
ggplot()+geom_point(data=C1, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#3E7C17')) + theme(legend.position = "none")

# C2
C2 <- subset(Data, Cluster=="C2")
ggplot()+geom_point(data=C2, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#3E7C27')) + theme(legend.position = "none")


