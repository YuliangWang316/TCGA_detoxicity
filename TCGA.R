Sys.setenv("VROOM_CONNECTION_SIZE" = 131072*5)
library(stringr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(dplyr)
library(org.Hs.eg.db)
library(tidyverse)
setwd("c:/Users/xjmik/Downloads/phaseI&phaseII_enzymes")
a<-list.files("c:/Users/xjmik/Downloads/phaseI&phaseII_enzymes")
d<-data.frame("Symbol")
colnames(d)<-"symbol"
for (i in 1:length(a)) {
  b<-read.csv(a[i],sep = ",",header = TRUE)
  c<-as.data.frame(b[,2])
  colnames(c)<-"symbol"
  d<-rbind(d,c)
}
remove(a,i,b,c)
d<-d[!duplicated(d$symbol),]
d<-d[2:506]

library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(EDASeq)
library(stringr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(org.Hs.eg.db)
library(tidyverse)
library(GSVA)
library(IDConverter)
library("limma")
setwd("D:/TEST3/")
TCGAtumor<-readRDS("d:/TCGAbiolinksPan_tumor.rds")
TCGAnormal<-readRDS("d:/TCGAbiolinksPan_normal.rds")
gene<-intersect(d,rownames(TCGAtumor))
TCGAnormal.data<-FindAllMarkers(TCGAnormal,only.pos = TRUE,logfc.threshold = 0,min.pct = 0,features = gene)
TCGAtumor.data<-FindAllMarkers(TCGAtumor,only.pos = TRUE,logfc.threshold = 0,min.pct = 0,features = gene)

TCGAnormal.data_filtered<-TCGAnormal.data[which(TCGAnormal.data$p_val_adj < 0.05),]
TCGAnormal.data_filtered$"PC.1-PCT.2"<-TCGAnormal.data_filtered$pct.1 - TCGAnormal.data_filtered$pct.2
TCGAnormal.data_filtered1<-TCGAnormal.data_filtered[which(TCGAnormal.data_filtered$`PC.1-PCT.2` > 0.5),]

TCGAnormal.data_filtered1 %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top10
a<-top10$gene
a<-unique(a)
DotPlot(TCGAnormal, features = a) 

TCGAtumor.data_filtered<-TCGAtumor.data[which(TCGAtumor.data$p_val_adj < 0.05),]
TCGAtumor.data_filtered$"PC.1-PCT.2"<-TCGAtumor.data_filtered$pct.1 - TCGAtumor.data_filtered$pct.2
TCGAtumor.data_filtered1<-TCGAtumor.data_filtered[which(TCGAtumor.data_filtered$`PC.1-PCT.2` > 0.5),]
TCGAtumor.data_filtered1 %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top10
a<-top10$gene
a<-unique(a)
DotPlot(TCGAtumor, features = a) 

DoHeatmap(TCGAnormal, features = gene ,size = 0,angle = 90) 
DoHeatmap(TCGAtumor, features = gene ,size = 0,angle = 90) 

Idents(TCGAnormal)<-TCGAnormal$group
Idents(TCGAtumor)<-TCGAtumor$group

LIHC_tumor<-subset(TCGAtumor,idents = "LIHC")
LIHC_normal<-subset(TCGAnormal,idents = "LIHC")

LIHC_count<-cbind(LIHC_tumor@assays$RNA@counts,LIHC_normal@assays$RNA@counts)  
LIHC_normal$GROUP<-rep("normal",length(colnames(LIHC_normal)))
LIHC_tumor$GROUP<-rep("tumor",length(colnames(LIHC_tumor)))

LIHC_metadata<-rbind(LIHC_tumor@meta.data[,c(1:5,40)],LIHC_normal@meta.data[,c(1:5,31)])  

pbmc <- CreateSeuratObject(counts = LIHC_count, project = "pbmc3k",meta.data = LIHC_metadata)  
remove(LIHC_count,LIHC_metadata)
remove(LIHC_normal,LIHC_tumor)

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)  
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)  
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))  
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)  
pbmc <- RunUMAP(pbmc, dims = 1:10)  
pbmc <- RunTSNE(pbmc, dims = 1:10)

Idents(pbmc)<-pbmc$GROUP
DimPlot(pbmc, reduction = "umap")  
DimPlot(pbmc, reduction = "tsne")
DoHeatmap(pbmc, features = gene ,size = 0,angle = 90)

gene<-intersect(d,rownames(pbmc))
remove(genename_1)
pbmc.tumor_vs_normal_detoxicity<-FindMarkers(pbmc,min.pct = 0,logfc.threshold = 0,ident.1 = "tumor",ident.2 = "normal",features = gene)  
remove(all.genes)

write.table(pbmc.markers,file = "c:/Users/xjmik/Desktop/LIHC_Allmarker_tumor_normal.txt",sep = "\t")
write.table(pbmc.markers_detoxicity,file = "c:/Users/xjmik/Desktop/LIHC_Allmarker_tumor_normal_detoxicity.txt",sep = "\t")
write.table(pbmc.tumor_vs_normal,file = "c:/Users/xjmik/Desktop/LIHC_tumor_vs_normal.txt",sep = "\t")
write.table(pbmc.tumor_vs_normal_detoxicity,file = "c:/Users/xjmik/Desktop/LIHC_tumor_vs_normal_detoxicity.txt",sep = "\t")
remove(pbmc.markers,pbmc.markers_detoxicity,pbmc.tumor_vs_normal,pbmc.tumor_vs_normal_detoxicity)

