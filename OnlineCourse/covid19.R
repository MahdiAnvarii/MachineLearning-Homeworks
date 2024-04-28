setwd("C:/Users/M/Desktop/Hi")
setRepositories()
library(data.table)
library(ggplot2)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2) # used for sequencing method
data <- fread("GSE147507_RawReadCounts_Ferret.tsv",data.table = F)
data[,-1] <- log2(data[,-1]+1)
data.m <- melt(data)
colnames(data.m) = c("Gene","Sample","Expression")
head(data)
rownames(data)=data$V1
data <- data[,-1]
ggplot(data,aes(x=Series10_FerretNW_Ctl_d1_1,y=Series10_FerretNW_Ctl_d1_2))+geom_point()
ggplot(data.m,aes(Sample,Expression))+geom_violin()
ggplot(data,aes(Series10_FerretNW_Ctl_d1_1))+geom_density()
ggplot(data.m,aes(Sample,Expression))+geom_boxplot(width=0.1,fill="White")

pdf("violine.pdf",width=15)
ggplot(data.m,aes(Sample,Expression,fill=Sample))+geom_violin()+geom_boxplot(width=0.1,fill="White")
dev.off()

col.sums <- colSums(data[,-1])
# normalizing with scaling for average
data.norm <- t(t(data[,-1])/col.sums)
colSums(data.norm)
colMeans(data.norm)
boxplot(data.norm)

# normalizing with DESeq2
data <- fread("GSE147507_RawReadCounts_Ferret.tsv",data.table = F)
gr <- c(rep("FerretNW_Ctl_d1",2),rep("FerretNW_SARS-CoV-2_d1",2),rep("FerretNW_Ctl_d3",2),rep("FerretNW_SARS-CoV-2_d3",2),rep("FerretNW_Ctl_d7",2),rep("FerretNW_SARS-CoV-2_d7",2),rep("FerretNW_IAV_d7",2),rep("FerretNW_Ctl_d14",2),rep("FerretNW_SARS-CoV-2_d14",2),rep("FerretTrachea_Ctl_d3",4),rep("FerretTrachea_IAV_d3",6),rep("FerretTrachea_SARS-CoV-2_d3",4))
gr
cds <- DESeqDataSetFromMatrix(data[,-1],colData = data.frame(Group=gr),~Group)
cds <- DESeq(cds)
data.norm <- log2(1+counts(cds,normalized=T))
head(data.norm)
boxplot(data.norm)

# dimensionality reduction by PCA
pc <- prcomp(data.norm)
plot(pc)
head(pc)
dim(pc$x)
head(pc$rotation)
dim(pc$rotation)
pcr <- data.frame(pc$r)
head(pcr)
dim(pcr)
pcr$group = gr
ggplot(pcr,aes(PC1,PC2,color=group))+geom_point(size=3 ,alpha=.8)+theme_bw()

pcx <- data.frame(pc$x)
head(pcx)
dim(pcx)
ggplot(pcx,aes(PC1,PC2))+geom_point()

# Destroy mean effect
data.mean.center <- t(scale(t(data.norm),scale=F))
head(data.mean.center)
head(rowMeans(data.mean.center))

pc2 <- prcomp(data.mean.center)
pcr2 <- data.frame(pc2$r)
head(pcr2)
dim(pcr2)
pcr2$group = gr
ggplot(pcr2,aes(PC1,PC2,color=group))+geom_point(size=3 ,alpha=.8)+theme_bw()

pcx2 <- data.frame(pc2$x)
head(pcx2)
dim(pcx2)
ggplot(pcx2,aes(PC1,PC2))+geom_point()
