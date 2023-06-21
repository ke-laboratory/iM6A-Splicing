rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(ggpubr)

#### Human
## Scatter plot
# Start
Data <- read.csv("iM6A_Human_IL_Start.csv")
Data <- Data[, c("chrom","name","strand","ID","m6AStart","m6AEnd","WT","IL","m6APosition")]
Data <- unique(Data)

Cluster <- read.csv("Human_Cluster_InternalExon.csv")
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
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
C1 <- subset(C1, m6APosition<=500)
table(C1$type)

# C2
C2 <- subset(Data, Cluster=="C2")
ggplot()+geom_point(data=C2, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
C2 <- subset(C2, m6APosition<=500)
table(C2$type)

### Scatter plot for <200nt
Data <- read.csv("iM6A_Human_IL_Start.csv")
Data <- Data[, c("chrom","name","strand","m6AStart","m6AEnd","WT","IL","m6APosition","ExonLength","ID")]
Data <- unique(Data)

Cluster <- read.csv("Human_Cluster_InternalExon.csv")
Data <- merge(Data, Cluster, by="ID")

Data <- subset(Data, ExonLength<200)

Gene <- Data[, c("chrom","name","strand")]
Gene <- unique(Gene)
Data <- merge(Data, Gene, by=c("chrom","name","strand"))

Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
table(Data$type)

Cluster <- Data[, c("ID","Cluster")]
C1 <- subset(Cluster, Cluster=="C1")
C2 <- subset(Cluster, Cluster=="C2")

C1 <- unique(C1)
C2 <- unique(C2)
set.seed(1)
C1 <- C1 %>% sample_n(3000, replace =FALSE)
set.seed(1)
C2 <- C2 %>% sample_n(3000, replace =FALSE)

Cluster <- rbind(C1,C2)

Data <- merge(Data, Cluster, by=c("ID","Cluster"))
Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"

# C1
C1 <- subset(Data, Cluster=="C1")
ggplot()+geom_point(data=C1, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.5) + ylim(-1,1) + xlim(0,200) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
table(C1$type)

# C2
C2 <- subset(Data, Cluster=="C2")
ggplot()+geom_point(data=C2, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.5) + ylim(-1,1) + xlim(0,200) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
table(C2$type)

### Scatter plot for 200-400nt
Data <- read.csv("iM6A_Human_IL_Start.csv")
Data <- Data[, c("chrom","name","strand","m6AStart","m6AEnd","WT","IL","m6APosition","ExonLength","ID")]
Data <- unique(Data)

Cluster <- read.csv("Human_Cluster_InternalExon.csv")
Data <- merge(Data, Cluster, by="ID")

Data <- subset(Data, ExonLength>=200&ExonLength<=400)

Gene <- Data[, c("chrom","name","strand")]
Gene <- unique(Gene)
Data <- merge(Data, Gene, by=c("chrom","name","strand"))

Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
table(Data$type)

Cluster <- Data[, c("ID","Cluster")]
C1 <- subset(Cluster, Cluster=="C1")
C2 <- subset(Cluster, Cluster=="C2")

C1 <- unique(C1)
C2 <- unique(C2)
set.seed(1)
C1 <- C1 %>% sample_n(1800, replace =FALSE)
set.seed(1)
C2 <- C2 %>% sample_n(1800, replace =FALSE)

Cluster <- rbind(C1,C2)

Data <- merge(Data, Cluster, by=c("ID","Cluster"))
Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"

# C1
C1 <- subset(Data, Cluster=="C1")
ggplot()+geom_point(data=C1, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.5) + ylim(-1,1) + xlim(0,400) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
table(C1$type)

# C2
C2 <- subset(Data, Cluster=="C2")
ggplot()+geom_point(data=C2, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.5) + ylim(-1,1) + xlim(0,400) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
table(C2$type)

### Scatter plot for >400nt
Data <- read.csv("iM6A_Human_IL_Start.csv")
Data <- Data[, c("chrom","name","strand","m6AStart","m6AEnd","WT","IL","m6APosition","ExonLength","ID")]
Data <- unique(Data)

Cluster <- read.csv("Human_Cluster_InternalExon.csv")
Data <- merge(Data, Cluster, by="ID")

Data <- subset(Data, ExonLength>400)

Gene <- Data[, c("chrom","name","strand")]
Gene <- unique(Gene)
Data <- merge(Data, Gene, by=c("chrom","name","strand"))

Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
table(Data$type)

Cluster <- Data[, c("ID","Cluster")]
C1 <- subset(Cluster, Cluster=="C1")
C2 <- subset(Cluster, Cluster=="C2")

C1 <- unique(C1)
C2 <- unique(C2)
set.seed(1)
C1 <- C1 %>% sample_n(150, replace =FALSE)
set.seed(1)
C2 <- C2 %>% sample_n(2000, replace =FALSE)

Cluster <- rbind(C1,C2)

Data <- merge(Data, Cluster, by=c("ID","Cluster"))
Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"

# C1
C1 <- subset(Data, Cluster=="C1")
ggplot()+geom_point(data=C1, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=1) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
C1 <- subset(C1, m6APosition<=500)
table(C1$type)


# C2
C2 <- subset(Data, Cluster=="C2")
ggplot()+geom_point(data=C2, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=1) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
C2 <- subset(C2, m6APosition<=500)
table(C2$type)


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
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
C1 <- subset(C1, m6APosition<=500)
table(C1$type)

# C2
C2 <- subset(Data, Cluster=="C2")
ggplot()+geom_point(data=C2, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.1) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
C2 <- subset(C2, m6APosition<=500)
table(C2$type)


### Scatter plot for <200nt
Data <- read.csv("iM6A_Mouse_IL_Start.csv")
Data <- Data[, c("chrom","name","strand","m6AStart","m6AEnd","WT","IL","m6APosition","ExonLength","ID")]
Data <- unique(Data)

Cluster <- read.csv("Mouse_Cluster_InternalExon.csv")
Data <- merge(Data, Cluster, by="ID")

Data <- subset(Data, ExonLength<200)

Gene <- Data[, c("chrom","name","strand")]
Gene <- unique(Gene)
Data <- merge(Data, Gene, by=c("chrom","name","strand"))

Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
table(Data$type)

Cluster <- Data[, c("ID","Cluster")]
C1 <- subset(Cluster, Cluster=="C1")
C2 <- subset(Cluster, Cluster=="C2")

C1 <- unique(C1)
C2 <- unique(C2)
set.seed(1)
C1 <- C1 %>% sample_n(3000, replace =FALSE)
set.seed(1)
C2 <- C2 %>% sample_n(3000, replace =FALSE)

Cluster <- rbind(C1,C2)

Data <- merge(Data, Cluster, by=c("ID","Cluster"))
Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"

# C1
C1 <- subset(Data, Cluster=="C1")
ggplot()+geom_point(data=C1, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.5) + ylim(-1,1) + xlim(0,200) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
table(C1$type)

# C2
C2 <- subset(Data, Cluster=="C2")
ggplot()+geom_point(data=C2, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.5) + ylim(-1,1) + xlim(0,200) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
table(C2$type)

### Scatter plot for 200-400nt
Data <- read.csv("iM6A_Mouse_IL_Start.csv")
Data <- Data[, c("chrom","name","strand","m6AStart","m6AEnd","WT","IL","m6APosition","ExonLength","ID")]
Data <- unique(Data)

Cluster <- read.csv("Mouse_Cluster_InternalExon.csv")
Data <- merge(Data, Cluster, by="ID")

Data <- subset(Data, ExonLength>=200&ExonLength<=400)

Gene <- Data[, c("chrom","name","strand")]
Gene <- unique(Gene)
Data <- merge(Data, Gene, by=c("chrom","name","strand"))

Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
table(Data$type)

Cluster <- Data[, c("ID","Cluster")]
C1 <- subset(Cluster, Cluster=="C1")
C2 <- subset(Cluster, Cluster=="C2")

C1 <- unique(C1)
C2 <- unique(C2)
set.seed(1)
C1 <- C1 %>% sample_n(1800, replace =FALSE)
set.seed(1)
C2 <- C2 %>% sample_n(1800, replace =FALSE)

Cluster <- rbind(C1,C2)

Data <- merge(Data, Cluster, by=c("ID","Cluster"))
Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"

# C1
C1 <- subset(Data, Cluster=="C1")
ggplot()+geom_point(data=C1, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.5) + ylim(-1,1) + xlim(0,400) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
table(C1$type)

# C2
C2 <- subset(Data, Cluster=="C2")
ggplot()+geom_point(data=C2, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=0.5) + ylim(-1,1) + xlim(0,400) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
table(C2$type)

### Scatter plot for >400nt
Data <- read.csv("iM6A_Mouse_IL_Start.csv")
Data <- Data[, c("chrom","name","strand","m6AStart","m6AEnd","WT","IL","m6APosition","ExonLength","ID")]
Data <- unique(Data)

Cluster <- read.csv("Mouse_Cluster_InternalExon.csv")
Data <- merge(Data, Cluster, by="ID")

Data <- subset(Data, ExonLength>400)

Gene <- Data[, c("chrom","name","strand")]
Gene <- unique(Gene)
Data <- merge(Data, Gene, by=c("chrom","name","strand"))

Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
table(Data$type)

Cluster <- Data[, c("ID","Cluster")]
C1 <- subset(Cluster, Cluster=="C1")
C2 <- subset(Cluster, Cluster=="C2")

C1 <- unique(C1)
C2 <- unique(C2)
set.seed(1)
C1 <- C1 %>% sample_n(165, replace =FALSE)
set.seed(1)
C2 <- C2 %>% sample_n(2000, replace =FALSE)

Cluster <- rbind(C1,C2)

Data <- merge(Data, Cluster, by=c("ID","Cluster"))
Data$Dvalue <- Data$IL - Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"

# C1
C1 <- subset(Data, Cluster=="C1")
ggplot()+geom_point(data=C1, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=1) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
C1 <- subset(C1, m6APosition<=500)
table(C1$type)

# C2
C2 <- subset(Data, Cluster=="C2")
ggplot()+geom_point(data=C2, aes(x=m6APosition, y=Dvalue, group=m6APosition, color=type), size=0.25, alpha=1) + ylim(-1,1) + xlim(0,500) +
  scale_color_manual(name = '', values = c('Up'='#DA1212','Non'='grey','Down'='#0000CD')) + theme(legend.position = "none")
C2 <- subset(C2, m6APosition<=500)
table(C2$type)






