library(plyr)
library(dplyr)
library(cluster) 
library(ggplot2)
library(pheatmap)
library(factoextra)

## Kmeans
df <- read.csv("mm10_Heatmap_LongInternalExonStart_200.csv", row.names = 1)
df = t(df)
df = data.frame(df)

fviz_nbclust(df, kmeans, method = "wss", k.max=10, iter.max = 100)
K <- kmeans(df, centers=2, nstart = 25, iter.max = 100)
fviz_cluster(object=K, data=df,
             ellipse.type = "euclid",star.plot=T,repel=T,
             geom = ("point"),palette='jco',main="",
             ggtheme=theme_minimal())+
  theme(axis.title = element_blank())


# Cluster info
Cluster <- data.frame(K$cluster)
Cluster$ID <- row.names(Cluster)
colnames(Cluster)[1] <- "Cluster"
Cluster <- Cluster[, c(2,1)]
Cluster <- Cluster[order(Cluster$Cluster, decreasing = FALSE),]

# Split data by Cluster
df$ID <- row.names(df)
df <- merge(df, Cluster, by=c("ID"))
df <- df[order(df$Cluster, decreasing = FALSE),]
row.names(df) <- df$ID

my_sample_col <- df[, c(1,42)]
my_sample_col$C[my_sample_col$Cluster==1] <- "C1"
my_sample_col$C[my_sample_col$Cluster==2] <- "C2"
my_sample_col <- my_sample_col[,-1]


df <- df[, -c(1,42)]
df <- t(df)

bk <- c(seq(-1,-0.001,by=0.001),seq(0,1,by=0.001))
p <- pheatmap(df, cluster_rows = FALSE, cluster_cols = FALSE,
              show_rownames = F, show_colnames = F,
              color = c(colorRampPalette(colors = c("green","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
              clustering_method = "complete",
              clustering_distance_cols = "euclidean",
              border_color = "white",
              cellwidth=0.05,
              cellheight=5,
              breaks=bk,
              annotation_col = my_sample_col)


# Exprot cluster
Cluster$Cluster[Cluster$Cluster==1] <- "C1"
Cluster$Cluster[Cluster$Cluster==2] <- "C2"
colnames(Cluster)[1] <- "ID"
write.csv(Cluster, "Mouse_Cluster_InternalExon.csv", row.names = F)


# Heatmap of Delta
Cluster <- read.csv("Mouse_Cluster_InternalExon.csv")
df <- read.csv("mm10_Heatmap_LongInternalExonStart_200.csv", row.names = 1)
df = t(df)
df = data.frame(df)
df$ID <- row.names(df)

df <- merge(df, Cluster, by="ID")
df <- df[order(df$Cluster, decreasing = FALSE),]
row.names(df) <- df$ID

my_sample_col <- df[, c(1,42)]
my_sample_col$C <- my_sample_col$Cluster
my_sample_col <- my_sample_col[,-1]

df <- df[, -c(1,42)]
df <- t(df)

bk <- c(seq(-1,-0.001,by=0.001),seq(0,1,by=0.001))
p <- pheatmap(df, cluster_rows = FALSE, cluster_cols = FALSE,
              show_rownames = F, show_colnames = F,
              color = c(colorRampPalette(colors = c("green","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
              clustering_method = "complete",
              clustering_distance_cols = "euclidean",
              border_color = "white",
              cellwidth=0.05,
              cellheight=6,
              breaks=bk,
              annotation_col = my_sample_col)


# Heatmap of iM6A
Cluster <- read.csv("Mouse_Cluster_InternalExon.csv")
df <- read.csv("mm10_Heatmap_iM6A_WT_LongInternalExonStart_200.csv", row.names = 1)
df = t(df)
df = data.frame(df)
df$ID <- row.names(df)

df <- merge(df, Cluster, by="ID")
df <- df[order(df$Cluster, decreasing = FALSE),]
row.names(df) <- df$ID

my_sample_col <- df[, c(1,42)]
my_sample_col$C <- my_sample_col$Cluster
my_sample_col <- my_sample_col[,-1]

df <- df[, -c(1,42)]
df <- t(df)

bk <- c(seq(0,0.4,by=0.001))
p <- pheatmap(df, cluster_rows = FALSE, cluster_cols = FALSE,
              show_rownames = F, show_colnames = F,
              color = c(colorRampPalette(colors = c("white","blue"))(length(bk)/2)),
              clustering_method = "complete",
              clustering_distance_cols = "euclidean",
              border_color = "white",
              cellwidth=0.05,
              cellheight=6,
              breaks=bk,
              annotation_col = my_sample_col)


Cluster <- read.csv("Mouse_Cluster_InternalExon.csv")
df <- read.csv("mm10_Heatmap_iM6A_IL_LongInternalExonStart_200.csv", row.names = 1)
df = t(df)
df = data.frame(df)
df$ID <- row.names(df)

df <- merge(df, Cluster, by="ID")
df <- df[order(df$Cluster, decreasing = FALSE),]
row.names(df) <- df$ID

my_sample_col <- df[, c(1,42)]
my_sample_col$C <- my_sample_col$Cluster
my_sample_col <- my_sample_col[,-1]

df <- df[, -c(1,42)]
df <- t(df)

bk <- c(seq(0,0.4,by=0.001))
p <- pheatmap(df, cluster_rows = FALSE, cluster_cols = FALSE,
              show_rownames = F, show_colnames = F,
              color = c(colorRampPalette(colors = c("white","blue"))(length(bk)/2)),
              clustering_method = "complete",
              clustering_distance_cols = "euclidean",
              border_color = "white",
              cellwidth=0.05,
              cellheight=6,
              breaks=bk,
              annotation_col = my_sample_col)


# Heatmap of RAC
Cluster <- read.csv("Mouse_Cluster_InternalExon.csv")
df <- read.csv("mm10_Heatmap_RAC_LongInternalExonStart_200.csv", row.names = 1)
df = t(df)
df = data.frame(df)
df$ID <- row.names(df)

df <- merge(df, Cluster, by="ID")
df <- df[order(df$Cluster, decreasing = FALSE),]
row.names(df) <- df$ID

my_sample_col <- df[, c(1,42)]
my_sample_col$C <- my_sample_col$Cluster
my_sample_col <- my_sample_col[,-1]

df <- df[, -c(1,42)]
df <- t(df)

bk <- c(seq(0,0.4,by=0.001))
p <- pheatmap(df, cluster_rows = FALSE, cluster_cols = FALSE,
              show_rownames = F, show_colnames = F,
              color = c(colorRampPalette(colors = c("white","blue"))(length(bk)/2)),
              clustering_method = "complete",
              clustering_distance_cols = "euclidean",
              border_color = "white",
              cellwidth=0.05,
              cellheight=6,
              breaks=bk,
              annotation_col = my_sample_col)


























