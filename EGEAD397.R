setwd("D:/E-GEAD-397.processed_autoimmunedisease_RNAseq/")
a = list.files("counts/")
dir = paste("./counts/",a,sep="")
  
library(dplyr)
library(tidyverse)
library(DESeq2)

b<-read.table(file = dir[1],sep="\t",row.names=1,header = TRUE)
b<-b[match(unique(b[,1]),b$Gene_name), ]
rownames(b)<-b[,1]
b<-b[,-1]
c<-as.data.frame(t(read.table("D:/E-GEAD-397.processed_autoimmunedisease_RNAseq/clinical_diagnosis_age_sex_v2.txt",sep = "\t",header = TRUE,row.names = 1)))
c<-c[,match(colnames(b),colnames(c))]
z <- strsplit(a[1],'-')
c[5,]<-rep(z[[1]][1],length(colnames(c)))
rownames(c)[5]<-"celltype"
d<-rbind(c,b)

for (j in 1:length(colnames(d))) {
  colnames(d)[j] <- paste(d[1,j],d[5,j],colnames(d)[j],sep = "-")
}
d<-d[-1,]
d<-d[-1,]
d<-d[-1,]
d<-d[-1,]
d<-d[-1,]
remove(j)
for (j in 1:length(colnames(c))) {
  colnames(c)[j] <- paste(c[1,j],c[5,j],colnames(c)[j],sep = "-")
}
c<-as.data.frame(t(c))
remove(b,z,j)
gene<-rownames(d)

for ( k in c(2:length(a))) {
  b<-read.table(file = dir[k],sep="\t",row.names=1,header = TRUE)
  b<-b[match(unique(b[,1]),b$Gene_name), ]
  rownames(b)<-b[,1]
  b<-b[,-1]
  e<-as.data.frame(t(read.table("D:/E-GEAD-397.processed_autoimmunedisease_RNAseq/clinical_diagnosis_age_sex_v2.txt",sep = "\t",header = TRUE,row.names = 1)))
  e<-e[,match(colnames(b),colnames(e))]
  z <- strsplit(a[k],'-')
  e[5,]<-rep(z[[1]][1],length(colnames(e)))
  rownames(e)[5]<-"celltype"
  f<-rbind(e,b)
  for (j in 1:length(colnames(f))) {
    colnames(f)[j] <- paste(f[1,j],f[5,j],colnames(f)[j],sep = "-")
  }
  remove(j)
  f<-f[-1,]
  f<-f[-1,]
  f<-f[-1,]
  f<-f[-1,]
  f<-f[-1,]
  for (j in 1:length(colnames(e))) {
    colnames(e)[j] <- paste(e[1,j],e[5,j],colnames(e)[j],sep = "-")
  }
  e<-as.data.frame(t(e))
  remove(b,z,j)
  gene2<-rownames(f)
  gene<-intersect(gene,gene2)
  remove(gene2)
  d<-d[gene,]
  f<-f[gene,]
  d<-cbind(d,f)
  remove(f)
  c<-rbind(c,e)
  remove(e)
}
remove(gene,k,a,dir)

library(dplyr)
library(Seurat)
library(patchwork)
pbmc <- CreateSeuratObject(counts = d, project = "pbmc3k",meta.data = c)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
remove(c,d,all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc,dims = 1:10)
pbmc <- RunLDA(pbmc,labels = pbmc$disease)
pbmc <- RunLDA(pbmc,labels = pbmc$celltype,reduction.name = "lda2")
pbmc <- RunUMAP(pbmc, dims = 1:10,reduction = "lda",reduction.name = "lda_umap")
pbmc <- RunTSNE(pbmc, dims = 1:10,reduction = "lda",reduction.name = "lda_tsne")
pbmc <- RunUMAP(pbmc, dims = 1:27,reduction = "lda2",reduction.name = "lda2_umap")
pbmc <- RunTSNE(pbmc, dims = 1:27,reduction = "lda2",reduction.name = "lda2_tsne")
Idents(pbmc)<-pbmc$celltype
DimPlot(pbmc,reduction = "lda2_umap")
DimPlot(pbmc,reduction = "lda2_tsne")
Idents(pbmc)<-pbmc$disease
DimPlot(pbmc,reduction = "lda_umap")
DimPlot(pbmc,reduction = "lda_tsne")
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
d<-intersect(d,rownames(pbmc))
Idents(pbmc)<-pbmc$celltype
pbmc_cell_markers<-FindAllMarkers(pbmc,only.pos = TRUE)
Idents(pbmc)<-pbmc$disease
pbmc_disease_markers<-FindAllMarkers(pbmc,only.pos = TRUE)
write.table(pbmc_disease_markers,file = "c:/Users/xjmik/Desktop/EGEAD397pbmc_disease_markers.txt",sep = "\t")
write.table(pbmc_cell_markers,file = "c:/Users/xjmik/Desktop/EGEAD397pbmc_cell_markers.txt",sep = "\t")
Idents(pbmc)<-pbmc$celltype
pbmc_cell_markers_detoxicity<-FindAllMarkers(pbmc,only.pos = TRUE,features = d,logfc.threshold = 0,min.pct = 0)
Idents(pbmc)<-pbmc$disease
pbmc_disease_markers_detoxicity<-FindAllMarkers(pbmc,only.pos = TRUE,features = d,logfc.threshold = 0,min.pct = 0)
saveRDS(pbmc,file = "c:/Users/xjmik/Desktop/EGEAD397_pbmc.rds")
write.table(pbmc_disease_markers_detoxicity,file = "c:/Users/xjmik/Desktop/EGEAD397pbmc_disease_markers_detoxicity.txt",sep = "\t")
write.table(pbmc_cell_markers_detoxicity,file = "c:/Users/xjmik/Desktop/EGEAD397pbmc_cell_markers_detoxicity.txt",sep = "\t")
remove(pbmc_cell_markers,pbmc_cell_markers_detoxicity,pbmc_disease_markers,pbmc_disease_markers_detoxicity)
Idents(pbmc)<-pbmc$celltype
pbmc_cell_markers_detoxicity %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
g1<-DoHeatmap(pbmc, features = top10$gene) 
Idents(pbmc)<-pbmc$disease
  pbmc_disease_markers_detoxicity %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
g2<-DoHeatmap(pbmc, features = top10$gene) 
ggsave(g2,filename = "c:/Users/xjmik/Desktop/Detoxicity_2/EGEAD397_disease.pdf",width = 10,height = 10)
ggsave(g1,filename = "c:/Users/xjmik/Desktop/Detoxicity_2/EGEAD397_celltype.pdf",width = 20,height = 20)
