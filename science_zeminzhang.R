PanT<-readRDS("D:/PanT.rds")
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
d<-intersect(d,rownames(PanT))
DimPlot(PanT)
PanT<-RunLDA(PanT,labels = PanT$meta.cluster)
PanT<-RunLDA(PanT,labels = PanT$cancerType,reduction.name = "lda2")
PanT<-RunLDA(PanT,labels = PanT$loc,reduction.name = "lda3")
PanT<-RunUMAP(PanT,reduction = "lda",reduction.name = "lda_umap",dims = 1:40)
PanT<-RunTSNE(PanT,reduction = "lda",reduction.name = "lda_tsne",dims = 1:40)
PanT<-RunUMAP(PanT,reduction = "lda2",reduction.name = "lda2_umap",dims = 1:17)
PanT<-RunTSNE(PanT,reduction = "lda2",reduction.name = "lda2_tsne",dims = 1:17)
PanT<-RunUMAP(PanT,reduction = "lda3",reduction.name = "lda3_umap",dims = 1:3)
PanT<-RunTSNE(PanT,reduction = "lda3",reduction.name = "lda3_tsne",dims = 1:3)
saveRDS(PanT,"d:/PanT.rds")