###########---------------------
rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(ggpubr)

Data <- read.csv("iM6A_Mouse_LastExon_Phylop.csv")
Data$Dvalue <- Data$Mutation-Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"
Data <- subset(Data, type != "Down")

# <=100
Data1 <- subset(Data, Data$m6APosition<=100)
x=colnames(Data1)[11]
y=colnames(Data1)[9]
group=levels(factor(Data1$type))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Data1, x="type", y="Score", fill = "type", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test")

Data1$type <- factor(Data1$type, levels=c("Up","Non"))
ggplot(Data1, aes(x=type, y=Score, fill=type)) + geom_boxplot(outlier.colour=NA) + ylim(-8,8)
table(Data1$type)
Up <- subset(Data1, type=="Up")
Non <- subset(Data1, type=="Non")
mean(Up$Score)
mean(Non$Score)

# >100
# Data1 <- subset(Data, Data$m6APosition>100 & Data$m6APosition<=200)
Data1 <- subset(Data, Data$m6APosition>100)
x=colnames(Data1)[11]
y=colnames(Data1)[9]
group=levels(factor(Data1$type))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Data1, x="type", y="Score", fill = "type", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test")

Data1$type <- factor(Data1$type, levels=c("Up","Non"))
ggplot(Data1, aes(x=type, y=Score, fill=type)) + geom_boxplot(outlier.colour=NA) + ylim(-8,8)
table(Data1$type)

Up <- subset(Data1, type=="Up")
Non <- subset(Data1, type=="Non")
mean(Up$Score)
mean(Non$Score)

###########---------------------
rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(ggpubr)

Data <- read.csv("iM6A_Mouse_SecondToLastExon_Phylop.csv")
Data$Dvalue <- Data$Deletion-Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"

iM6A <- read.csv("iM6A_Mouse_SecondToLastExon.csv")
iM6A <- iM6A[, c(1,4,5,9)]
colnames(iM6A) <- c("chrom","name","strand","ExonLength")
iM6A <- unique(iM6A)

Data <- merge(Data, iM6A, by=c("chrom","name","strand"))
Data$m6APosition <- Data$ExonLength - Data$Distance

Data <- subset(Data, type != "Down")

# <=100
Data1 <- subset(Data, Data$m6APosition<=100)
x=colnames(Data1)[11]
y=colnames(Data1)[9]
group=levels(factor(Data1$type))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Data1, x="type", y="Score", fill = "type", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test")

Data1$type <- factor(Data1$type, levels=c("Up","Non"))
ggplot(Data1, aes(x=type, y=Score, fill=type)) + geom_boxplot(outlier.colour=NA) + ylim(-8,8)
table(Data1$type)

Up <- subset(Data1, type=="Up")
Non <- subset(Data1, type=="Non")
mean(Up$Score)
mean(Non$Score)

# >100
Data1 <- subset(Data, Data$m6APosition>100 & Data$Distance>100 & Data$m6APosition<=200)
x=colnames(Data1)[11]
y=colnames(Data1)[9]
group=levels(factor(Data1$type))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Data1, x="type", y="Score", fill = "type", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test")

Data1$type <- factor(Data1$type, levels=c("Up","Non"))
ggplot(Data1, aes(x=type, y=Score, fill=type)) + geom_boxplot(outlier.colour=NA) + ylim(-8,8)
table(Data1$type)

Up <- subset(Data1, type=="Up")
Non <- subset(Data1, type=="Non")
mean(Up$Score)
mean(Non$Score)





###########---------------------
rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(ggpubr)

Data <- read.csv("iM6A_Mouse_InternalExon_Phylop.csv")
Data$Dvalue <- Data$IL-Data$WT
Data$type[Data$Dvalue > 0.1] <-"Up"
Data$type[Data$Dvalue < -0.1] <-"Down"
Data$type[Data$Dvalue <= 0.1 & Data$Dvalue >= -0.1] <-"Non"

iM6A <- read.csv("iM6A_Mouse_IL_Start.csv")
iM6A <- iM6A[, c("ID", "ExonLength")]
iM6A <- unique(iM6A)

Data <- merge(Data, iM6A, by=c("ID"))
Data$m6APosition1 <- Data$m6APosition
Data$m6APosition2 <- Data$ExonLength - Data$m6APosition

Data <- subset(Data, type != "Down")

# <=100
Data1 <- subset(Data, m6APosition1<=100 | m6APosition2<=100)
x=colnames(Data1)[11]
y=colnames(Data1)[9]
group=levels(factor(Data1$type))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Data1, x="type", y="Score", fill = "type", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test")

Data1$type <- factor(Data1$type, levels=c("Up","Non"))
ggplot(Data1, aes(x=type, y=Score, fill=type)) + geom_boxplot(outlier.colour=NA) + ylim(-8,8)
table(Data1$type)

Up <- subset(Data1, type=="Up")
Non <- subset(Data1, type=="Non")
mean(Up$Score)
mean(Non$Score)

# >100
Data1 <- subset(Data, m6APosition1>100 & m6APosition2>100)
x=colnames(Data1)[11]
y=colnames(Data1)[9]
group=levels(factor(Data1$type))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
ggboxplot(Data1, x="type", y="Score", fill = "type", 
          xlab=x, ylab=y,
          legend.title=x,
          add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test")

Data1$type <- factor(Data1$type, levels=c("Up","Non"))
ggplot(Data1, aes(x=type, y=Score, fill=type)) + geom_boxplot(outlier.colour=NA) + ylim(-8,8)
table(Data1$type)

Up <- subset(Data1, type=="Up")
Non <- subset(Data1, type=="Non")
mean(Up$Score)
mean(Non$Score)


