library(dplyr)
library(patchwork)
library(Seurat)

pbmc.data<-Read10X("c:/Users/xjmik/Downloads/GSE189175/DATA/")
pbmc.metadata<-read.table("c:/Users/xjmik/Downloads/GSE189175/meta_data.tsv",sep = "\t",header = TRUE,row.names = 1)

pbmc.data<-pbmc.data[na.omit(rownames(pbmc.data)),]
pbmc.metadata<-pbmc.metadata[colnames(pbmc.data),]
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",meta.data = pbmc.metadata)
remove(pbmc.data,pbmc.metadata)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
remove(all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 3)
pbmc <- RunUMAP(pbmc,dims = 1:30)
DimPlot(pbmc)
pbmc <- RunTSNE(pbmc,dims = 1:30)
DimPlot(pbmc,reduction = "tsne")
