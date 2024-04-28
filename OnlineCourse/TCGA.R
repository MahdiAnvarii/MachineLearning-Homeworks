# Clustering
setRepositories()
library(DESeq2)
library(cli)
library(TCGAbiolinks)
remove.packages("cli")
install.packages("cli")
install.packages("TCGAbiolinks")
library(tidyverse)
library(pheatmap)
library(maftools)
library(SummarizedExperiment)

gdcprojects <- getGDCprojects()
gdcprojects
getProjectSummary('TCGA-BRCA')

query_TCGA <- GDCquery(project = 'TCGA-BRCA' , data.category = 'Transcriptome Profiling')
output_query_TCGA <- getResults(query_TCGA)

query_TCGA <- GDCquery(project = 'TCGA-BRCA' , 
                       data.category = 'Transcriptome Profiling' ,
                       experimental.strategy = 'RNA-Seq' ,
                       access = 'open' ,
                       barcode = c('TCGA-BH-A1FU-11A-23R-A14D-07','TCGA-BH-A18M-11A-33R-A12D-07',
                                   'TCGA-BH-A18R-11A-42R-A12D-07','TCGA-BH-A1FE-11B-14R-A13Q-07',
                                   'TCGA-BH-A0BZ-11A-61R-A12P-07','TCGA-BH-A0DH-11A-31R-A089-07',
                                   'TCGA-BH-A1F0-11B-23R-A137-07','TCGA-BH-A0BS-11A-11R-A12P-07',
                                   'TCGA-A7-A26E-01B-06R-A277-07','TCGA-A2-A0CU-01A-12R-A034-07',
                                   'TCGA-PL-A8LV-01A-21R-A41B-07','TCGA-BH-A0BC-01A-22R-A084-07',
                                   'TCGA-AR-A1AX-01A-11R-A12P-07','TCGA-AC-A2FO-01A-11R-A180-07',
                                   'TCGA-AQ-A0Y5-01A-11R-A14M-07','TCGA-AC-A3EH-01A-22R-A22K-07'))
output_query_TCGA2 <- getResults(query_TCGA)

GDCdownload(query_TCGA)

tcga_brca_data <- GDCprepare(query_TCGA,summarizedExperiment = TRUE)
brca_matrix <- assay(tcga_brca_data,'fpkm_unstrand')

saveRDS(brca_matrix,"TCGA.rds")
getwd()

my_data <- readRDS("TCGA.rds")
my_data <- my_data[,sample(1:ncol(my_data))]

cl <- colnames(my_data)
cl <- as.array(cl)
jj <- sapply(1:ncol(my_data), function(i) str_detect(cl[i],'11A') | str_detect(cl[i],'11B'))
jj
tumor = cl[jj]
normal = cl[!jj]

gr=rep("X",16)
gr[jj]<- "T"
gr[!jj]<- "N"
gr
colnames(my_data) <- make.names(gr,unique = T)

data <- my_data[which(rowSums(my_data>2)>=5),]
dim(my_data)
dim(data)

boxplot(data)
boxplot(log2(1+data))

cold = data.frame(Group=gr)
rownames(cold) = colnames(data)
cds <- DESeqDataSetFromMatrix(round(data),colData = cold,~Group)
cds <- DESeq(cds)
data.norm <- log2(1+counts(cds,normalized=T))
boxplot(data.norm)

# hierarchical clustering
pheatmap(data.norm)

# mean center
pheatmap(t(scale(t(data.norm))))

pheatmap(cor(data.norm))

install.packages("gplots")
library(gplots)

annotation_for_heatmap <- data.frame(Disease = gr)
row.names(annotation_for_heatmap) <- colnames(data.norm)
pheatmap(cor(data.norm),labels_row = colnames(data.norm),
         annotation_row = annotation_for_heatmap, 
         labels_col =colnames(data.norm), color = bluered(256),border_color = NA ,
         main = "heatmap between samples")

# k means
vari <- apply(data.norm, 1, var)
hist(vari)
data.norm2 <- data.norm[vari>0.1,]

data.norm2 <- data.norm2-rowMeans(data.norm2)
clu <- kmeans(data.norm2,100)
clu
table(clu$cluster)

library(reshape)
exmf <- data.frame(Gene=rownames(data.norm2),data.norm2,Cluster=factor(clu$cluster))
mex <- reshape::melt(exmf)
head(mex)
colnames(mex)[3:4] <- c('Sample','Expression')

pdf("ResultClustering.pdf",width = 15, height = 10)
lapply(1:100,function(i){ggplot(subset(mex,Cluster==i),aes(Sample,Expression,group=Gene))+
    geom_line(alpha=0.4) + labs(x=paste0("Cluster:",i)) })
dev.off()

ggplot(subset(mex,Cluster==4),aes(Sample,Expression,group=Gene))+
  geom_line(alpha=0.4) + labs(x=paste0("Cluster:",4))