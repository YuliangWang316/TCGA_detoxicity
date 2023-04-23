library(dplyr)
library(Seurat)

pbmc.data<-read.table("c:/Users/xjmik/Downloads/GSE108989/GSE108989_CRC.TCell.S11138.count.txt",sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
pbmc.data<-pbmc.data[!duplicated(pbmc.data$symbol),]
pbmc.data<-na.omit(pbmc.data)
rownames(pbmc.data)<-pbmc.data$symbol
pbmc.data<-pbmc.data[,-1]
pbmc.metadata<-read.table("c:/Users/xjmik/Downloads/GSE108989/metadata.txt",sep = "\t",header = TRUE,row.names = 1)
pbmc.metadata<-pbmc.metadata[colnames(pbmc.data),]
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",meta.data = pbmc.metadata)
remove(pbmc.data,pbmc.metadata)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
remove(all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(pbmc,dims = 1:20)
Idents(pbmc)<-pbmc@meta.data$majorCluster

pbmc.marker<-FindAllMarkers(pbmc,only.pos = TRUE)
write.table(pbmc.marker,"C:/Users/xjmik/Desktop/GSE108989_Totalmarker.txt",sep = "\t")
remove(pbmc.marker)
for (i in 1:length(rownames(pbmc@meta.data))) {
  if(pbmc@meta.data$sampleType[i] == "PTH"){
    pbmc@meta.data$Bulksample[i] <- "P"
  }
  if(pbmc@meta.data$sampleType[i] == "PTC"){
    pbmc@meta.data$Bulksample[i] <- "P"
  }
  if(pbmc@meta.data$sampleType[i] == "PTR"){
    pbmc@meta.data$Bulksample[i] <- "P"
  }
  if(pbmc@meta.data$sampleType[i] == "TTH"){
    pbmc@meta.data$Bulksample[i] <- "T"
  }
  if(pbmc@meta.data$sampleType[i] == "TTR"){
    pbmc@meta.data$Bulksample[i] <- "T"
  }
  if(pbmc@meta.data$sampleType[i] == "TTC"){
    pbmc@meta.data$Bulksample[i] <- "T"
  }
  if(pbmc@meta.data$sampleType[i] == "NTY"){
    pbmc@meta.data$Bulksample[i] <- "N"
  }
  if(pbmc@meta.data$sampleType[i] == "TTY"){
    pbmc@meta.data$Bulksample[i] <- "T"
  }
  if(pbmc@meta.data$sampleType[i] == "PTY"){
    pbmc@meta.data$Bulksample[i] <- "P"
  }
  if(pbmc@meta.data$sampleType[i] == "NTH"){
    pbmc@meta.data$Bulksample[i] <- "N"
  }
  if(pbmc@meta.data$sampleType[i] == "NTR"){
    pbmc@meta.data$Bulksample[i] <- "N"
  }
  if(pbmc@meta.data$sampleType[i] == "NTC"){
    pbmc@meta.data$Bulksample[i] <- "N"
  }
  if(pbmc@meta.data$sampleType[i] == "NP7"){
    pbmc@meta.data$Bulksample[i] <- "N"
  }
  if(pbmc@meta.data$sampleType[i] == "PP7"){
    pbmc@meta.data$Bulksample[i] <- "P"
  }
  if(pbmc@meta.data$sampleType[i] == "TP7"){
    pbmc@meta.data$Bulksample[i] <- "T"
  }
}
remove(i)
for (i in levels(Idents(pbmc))[c(1:18,22:26,28)]) {
  pbmcs<-subset(pbmc,idents = i)
  Idents(pbmcs)<-pbmcs@meta.data$Bulksample
  pbmcs.markers<-FindMarkers(pbmcs,only.pos = TRUE,ident.1 = "T",ident.2 = "N")
  write.table(pbmcs.markers,file = paste0("c:/Users/xjmik/Desktop/GSE108989_T_N",i,".txt"),sep = "\t")
}
remove(i,pbmcs,pbmcs.markers)
for (i in levels(Idents(pbmc))[c(2:5,8,10,12:15,17,20:22,28)]) {
  pbmcs<-subset(pbmc,idents = i)
  Idents(pbmcs)<-pbmcs@meta.data$Bulksample
  pbmcs.markers<-FindMarkers(pbmcs,only.pos = TRUE,ident.1 = "T",ident.2 = "P")
  write.table(pbmcs.markers,file = paste0("c:/Users/xjmik/Desktop/GSE108989_T_P",i,".txt"),sep = "\t")
  remove(pbmcs.markers)
}

remove(i,pbmcs)
for (i in 1:length(rownames(pbmc@meta.data))) {
  if(pbmc@meta.data$sampleType[i] == "PTH"){
    pbmc@meta.data$Bulksample2[i] <- "TH"
  }
  if(pbmc@meta.data$sampleType[i] == "PTC"){
    pbmc@meta.data$Bulksample2[i] <- "TC"
  }
  if(pbmc@meta.data$sampleType[i] == "PTR"){
    pbmc@meta.data$Bulksample2[i] <- "TR"
  }
  if(pbmc@meta.data$sampleType[i] == "TTH"){
    pbmc@meta.data$Bulksample2[i] <- "TH"
  }
  if(pbmc@meta.data$sampleType[i] == "TTR"){
    pbmc@meta.data$Bulksample2[i] <- "TR"
  }
  if(pbmc@meta.data$sampleType[i] == "TTC"){
    pbmc@meta.data$Bulksample2[i] <- "TC"
  }
  if(pbmc@meta.data$sampleType[i] == "NTY"){
    pbmc@meta.data$Bulksample2[i] <- "TY"
  }
  if(pbmc@meta.data$sampleType[i] == "PTY"){
    pbmc@meta.data$Bulksample2[i] <- "TY"
  }
  if(pbmc@meta.data$sampleType[i] == "TTY"){
    pbmc@meta.data$Bulksample2[i] <- "TY"
  }
  if(pbmc@meta.data$sampleType[i] == "NTH"){
    pbmc@meta.data$Bulksample2[i] <- "TH"
  }
  if(pbmc@meta.data$sampleType[i] == "NTR"){
    pbmc@meta.data$Bulksample2[i] <- "TR"
  }
  if(pbmc@meta.data$sampleType[i] == "NTC"){
    pbmc@meta.data$Bulksample2[i] <- "TC"
  }
  if(pbmc@meta.data$sampleType[i] == "NP7"){
    pbmc@meta.data$Bulksample2[i] <- "P7"
  }
  if(pbmc@meta.data$sampleType[i] == "PP7"){
    pbmc@meta.data$Bulksample2[i] <- "P7"
  }
  if(pbmc@meta.data$sampleType[i] == "TP7"){
    pbmc@meta.data$Bulksample2[i] <- "P7"
  }
}
remove(i)
Idents(pbmc)<-pbmc@meta.data$Bulksample2
pbmc.marker2<-FindAllMarkers(pbmc,only.pos = TRUE)
write.table(pbmc.marker2,"C:/Users/xjmik/Desktop/GSE108989_BulkTotalmarkers.txt",sep = "\t")
remove(pbmc.marker2)
for (i in levels(Idents(pbmc))) {
  pbmcs<-subset(pbmc,idents = i)
  Idents(pbmcs)<-pbmcs@meta.data$Bulksample
  pbmcs.markers<-FindMarkers(pbmcs,only.pos = TRUE,ident.1 = "T",ident.2 = "N")
  write.table(pbmcs.markers,file = paste0("c:/Users/xjmik/Desktop/GSE108989_Bulk_T_N",i,".txt"),sep = "\t")
}
remove(i,pbmcs,pbmcs.markers)
for (i in levels(Idents(pbmc))) {
  pbmcs<-subset(pbmc,idents = i)
  Idents(pbmcs)<-pbmcs@meta.data$Bulksample
  pbmcs.markers<-FindMarkers(pbmcs,only.pos = TRUE,ident.1 = "T",ident.2 = "P")
  write.table(pbmcs.markers,file = paste0("c:/Users/xjmik/Desktop/GSE108989_Bulk_T_P",i,".txt"),sep = "\t")
  remove(pbmcs.markers)
}
remove(i,pbmcs)
pbmc<-RunLDA(pbmc,labels = pbmc$majorCluster)
pbmc<-RunLDA(pbmc,labels = pbmc$Bulksample2,reduction.name = "lda_bulksample2")
pbmc<-RunTSNE(pbmc,dims = 1:20)
pbmc<-RunICA(pbmc)
pbmc<-RunSLSI(pbmc,graph = "RNA_snn")
pbmc<-RunSPCA(pbmc,graph = "RNA_snn")
pbmc<-RunUMAP(pbmc,reduction = "lda",reduction.name = "lda_umap",dims = 1:28)
pbmc<-RunTSNE(pbmc,reduction = "lda",reduction.name = "lda_tsne",dims = 1:28)
Idents(pbmc)<-pbmc@meta.data$majorCluster
DimPlot(pbmc,reduction = "lda_umap")
DimPlot(pbmc,reduction = "lda_tsne")
pbmc<-RunUMAP(pbmc,reduction = "lda_bulksample2",reduction.name = "lda_bulksample2_umap",dims = 1:4)
pbmc<-RunTSNE(pbmc,reduction = "lda_bulksample2",reduction.name = "lda_bulksample2_tsne",dims = 1:4)
Idents(pbmc)<-pbmc@meta.data$Bulksample2
DimPlot(pbmc,reduction = "lda_bulksample2_umap")
DimPlot(pbmc,reduction = "lda_bulksample2_tsne")
saveRDS(pbmc,file = "c:/Users/xjmik/Desktop/GSE108989pbmc.rds")
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
Idents(pbmc)<-pbmc$majorCluster
pbmc.markers_detoxicity<-FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0,min.pct = 0,features = d)
write.table(pbmc.markers_detoxicity,file = "c:/Users/xjmik/Desktop/GSE108989marker_Detoxicity.txt",sep = "\t")
Idents(pbmc)<-pbmc$Bulksample
pbmc.markers_detoxicity_big<-FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0,min.pct = 0,features = d)
write.table(pbmc.markers_detoxicity_big,file = "c:/Users/xjmik/Desktop/GSE108989marker_Detoxicity_tissue.txt",sep = "\t")
Idents(pbmc)<-pbmc$majorCluster
pbmc.markers_detoxicity %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
g1<-DoHeatmap(pbmc, features = top10$gene)
library(ggplot2)
ggsave(g1,filename = "c:/Users/xjmik/Desktop/Detoxicity_2/GSE108989_celltype.pdf",width = 30,height = 20)
Idents(pbmc)<-pbmc$Bulksample
pbmc.markers_detoxicity_big %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
g2<-DoHeatmap(pbmc, features = top10$gene)
library(ggplot2)
ggsave(g2,filename = "c:/Users/xjmik/Desktop/Detoxicity_2/GSE108989_tissue.pdf",width = 10,height = 15)
Idents(pbmc)<-pbmc$majorCluster
FeaturePlot(pbmc,features = "HTATIP2",reduction = "lda_tsne")
