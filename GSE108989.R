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
DimPlot(pbmc,reduction = "lda_tsne")
Idents(pbmc)<-pbmc@meta.data$majorCluster
pbmc_avg<-AverageExpression(pbmc)
write.table(pbmc_avg$RNA,file = "c:/Users/xjmik/Downloads/GSE98638/average_Expression_cluster.txt",sep = "\t")
remove(pbmc_avg)
pbmc.markers<-FindMarkers(pbmc,only.pos = TRUE,ident.1 = "C08_CD4-CTLA4",ident.2 = "C07_CD4-FOXP3")
write.table(pbmc.markers,file = "c:/Users/xjmik/Desktop/GSE98638/Total_CTLA4_vs_FOXP3.txt",sep = "\t")
remove(pbmc.markers)
Idents(pbmc)<-pbmc@meta.data$Bulksample
DimPlot(pbmc,reduction = "lda_tsne")
FeaturePlot(pbmc,features = "CTLA4",reduction = "lda_tsne")
Idents(pbmc)<-pbmc@meta.data$majorCluster
DotPlot(pbmc, features = c("S100A4","DYNLL1","GPX1","ATP5E","MYL6")) + RotatedAxis()
venn<-read.table("c:/Users/xjmik/Desktop/GSE98638/Bulk3_venn_dotplot.txt",sep = "\t")
venn<-venn[,1]
DotPlot(pbmc, features = venn) + RotatedAxis()
remove(venn)
Idents(pbmc)<-pbmc$Bulksample
pbmcn<-subset(pbmc,idents = c("T","N","P"))
DotPlot(pbmcn, features = c("S100A4","DYNLL1","GPX1","ATP5E","MYL6")) + RotatedAxis()
venn<-read.table("c:/Users/xjmik/Desktop/GSE98638/Bulk3_venn_dotplot.txt",sep = "\t")
venn<-venn[,1]
DotPlot(pbmcn, features = venn) + RotatedAxis()
remove(venn,pbmcn)
Idents(pbmc)<-pbmc$majorCluster

h<-read.table("c:/Users/xjmik/Downloads/hallmarks.txt",sep = "\t")
rownames(h)<-h$V1
h<-as.data.frame(t(h))
h<-h[-1,]
h<-h[-1,]
h<-as.list(h)
library(GSVA)
gsva<-gsva(expr = as.matrix(pbmc@assays$RNA@data),gset.idx.list = h,kcdf="Poisson",parallel.sz=20,method ="ssgsea")
write.table(gsva,file="c:/Users/xjmik/Desktop/GSE98638/HallmarksGSVA.txt",sep="\t")

h<-read.table("c:/Users/xjmik/Downloads/c2CGP.txt",sep = "\t")
rownames(h)<-h$V1
h<-as.data.frame(t(h))
h<-h[-1,]
h<-h[-1,]
h<-as.list(h)
library(GSVA)
gsva_2<-gsva(expr = as.matrix(pbmc@assays$RNA@data),gset.idx.list = h,kcdf="Poisson",parallel.sz=20,method ="ssgsea")
write.table(gsva,file="c:/Users/xjmik/Desktop/GSE98638/c2CGPGSVA.txt",sep="\t")

h<-read.table("c:/Users/xjmik/Downloads/c2KEGG.txt",sep = "\t")
rownames(h)<-h$V1
h<-as.data.frame(t(h))
h<-h[-1,]
h<-h[-1,]
h<-as.list(h)
library(GSVA)
gsva_3<-gsva(expr = as.matrix(pbmc@assays$RNA@data),gset.idx.list = h,kcdf="Poisson",parallel.sz=20,method ="ssgsea")
write.table(gsva,file="c:/Users/xjmik/Desktop/GSE98638/c2KEGGGSVA.txt",sep="\t")

Idents(pbmc)<-pbmc@meta.data$majorCluster
Tregtumor<-subset(pbmc,idents = c("CD4_C12-CTLA4"))
Tregpbmc<-subset(pbmc,idents = c("CD4_C10-FOXP3"))
Tregpbmc.data<-Tregpbmc@assays$RNA@counts
remove(Tregpbmc)
JMJD1C<-FetchData(Tregtumor,vars = "JMJD1C")
JMJD1C$sample<-rep("jmjd1c",length(rownames(JMJD1C)))

Namejmjd1chi<-rownames(JMJD1C[which(JMJD1C$JMJD1C>mean(JMJD1C$JMJD1C)),])
Namejmjd1clo<-rownames(JMJD1C[which(JMJD1C$JMJD1C<=mean(JMJD1C$JMJD1C)),])
JMJD1Chi.data<-as.data.frame(Tregtumor@assays$RNA@counts)[,Namejmjd1chi]
JMJD1Clo.data<-as.data.frame(Tregtumor@assays$RNA@counts)[,Namejmjd1clo]
remove(Namejmjd1chi,Namejmjd1clo,Tregtumor,JMJD1C)

JMJD1Chi.data <- as.data.frame(JMJD1Chi.data)
JMJD1Clo.data <- as.data.frame(JMJD1Clo.data)
Tregpbmc.data <- as.data.frame(Tregpbmc.data)

for (i in 1:length(colnames(JMJD1Chi.data))) {
  colnames(JMJD1Chi.data)[i] <- paste(colnames(JMJD1Chi.data)[i],"JMJD1Chi",i,sep = "-")  
}

for (i in 1:length(colnames(JMJD1Clo.data))) {
  colnames(JMJD1Clo.data)[i] <- paste(colnames(JMJD1Clo.data)[i],"JMJD1Clo",i,sep = "-")  
}

for (i in 1:length(colnames(Tregpbmc.data))) {
  colnames(Tregpbmc.data)[i] <- paste(colnames(Tregpbmc.data)[i],"Tregpbmc",i,sep = "-")  
}

JMJD1Chi.metadata<-data.frame(colnames(JMJD1Chi.data),rep("hi",length(colnames(JMJD1Chi.data))))
JMJD1Clo.metadata<-data.frame(colnames(JMJD1Clo.data),rep("lo",length(colnames(JMJD1Clo.data))))
Tregpbmc.metadata<-data.frame(colnames(Tregpbmc.data),rep("pbmc",length(colnames(Tregpbmc.data))))
colnames(JMJD1Chi.metadata)<-c("barcode","group")
colnames(JMJD1Clo.metadata)<-c("barcode","group")
colnames(Tregpbmc.metadata)<-c("barcode","group")
pbmc.metadata<-rbind(JMJD1Chi.metadata,JMJD1Clo.metadata,Tregpbmc.metadata)
rownames(pbmc.metadata)<-pbmc.metadata[,1]
pbmc.data<-cbind(JMJD1Chi.data,JMJD1Clo.data,Tregpbmc.data)
remove(i,JMJD1Chi.metadata,JMJD1Clo.metadata,Tregpbmc.metadata,JMJD1Chi.data,JMJD1Clo.data,Tregpbmc.data)
pbmc_raw<-pbmc
remove(pbmc)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",meta.data = pbmc.metadata)
remove(pbmc.data,pbmc.metadata)

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(pbmc, dims = 1:20)
pbmc <- RunTSNE(pbmc, dims = 1:20)
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "umap",split.by = "group")
DimPlot(pbmc, reduction = "tsne")
DimPlot(pbmc, reduction = "tsne",split.by = "group")


Idents(pbmc)<-pbmc@meta.data$group
pbmc<-RunLDA(pbmc,labels = pbmc$group)
DimPlot(pbmc, reduction = "lda")
pbmc<-RunUMAP(pbmc,reduction = "lda",reduction.name = "lda_umap",dims = 1:2)
pbmc<-RunTSNE(pbmc,reduction = "lda",reduction.name = "lda_tsne",dims = 1:2)

DimPlot(pbmc,reduction = "lda_umap")
DimPlot(pbmc,reduction = "lda_tsne")
DimPlot(pbmc,reduction = "lda_umap")
pbmc_new_3<-as.SingleCellExperiment(pbmc)

library(slingshot)
sce_6<-slingshot(pbmc_new_3,clusterLabels = "group",reducedDim = 'LDA',start.clus = "pbmc",end.clus = "hi",dist.method= "slingshot") 
library(TrajectoryUtils)
library(grDevices)
library(RColorBrewer)
library(ggplot2)
library(ggsn)
library(scales)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce_6$slingPseudotime_1, breaks=100)]
plot(reducedDims(sce_6)$LDA, col = plotcol, pch=16, asp = 1) 
lines(SlingshotDataSet(sce_6), lwd=2, col='black')
data<-sce_6@int_colData$reducedDims@listData$LDA
color<-data.frame(sce_6$slingPseudotime_1,plotcol)
rownames(color)<-rownames(data)
data<-cbind(data,color)
colnames(data)[3]<-"Pseudotime"
colnames(data)[4]<-"color"
remove(color)
library(ggplot2)
library(scales)
library(ggsn)
ggplot(data, aes(x = LDA_1, y = LDA_2,fill=Pseudotime)) +
  geom_point(col = plotcol ) + 
  theme_classic() +
  scale_fill_gradientn(colors = colors)
remove(plotcol,colors,all.genes,sce_6,pbmc_new_3,data)

library(monocle)
trace('project2MST',edit = T,where = asNamespace("monocle"))
data<-as.sparse(pbmc@assays$RNA@counts)
pd <-pbmc@meta.data
pd <- new('AnnotatedDataFrame', data = pd)  
fData<-data.frame(gene_short_name = row.names(data),row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data, phenoData = pd,featureData = fd,expressionFamily=negbinomial())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
# 
# 
# 
pbmc.marker<-FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0,min.pct = 0)
# diff_test_res <- differentialGeneTest(monocle_cds,fullModelFormulaStr = "~group",cores = 20)
# write.table(diff_test_res,"c:/Users/xjmik/Desktop/difftest.txt",sep = "\t")
write.table(pbmc.marker,"c:/Users/xjmik/Desktop/pbmcmarkers_3.txt",sep = "\t")
pbmcmarkers<-read.table("c:/Users/xjmik/Desktop/pbmcmarkers_3.txt",sep = "\t",header = TRUE,row.names = 1)
# diff_test_res<-read.table("c:/Users/xjmik/Desktop/difftest.txt",sep = "\t",header = TRUE,row.names = 1)
# remove(diff_test_res)
pbmcmarkers_new<-pbmcmarkers[which(pbmcmarkers$p_val_adj < 0.05 & pbmcmarkers$avg_log2FC > 0.25 & pbmcmarkers$pct.1 > 0.1 & pbmcmarkers$pct.2 > 0.1),]
pbmcmarkers_new$filter<-pbmcmarkers_new$pct.1 - pbmcmarkers_new$pct.2
pbmcmarkers_new_new<-pbmcmarkers_new[which(pbmcmarkers_new$filter > median(pbmcmarkers_new$filter)),]
ordering_genes <- pbmcmarkers_new_new$gene
monocle_cds <-setOrderingFilter(monocle_cds,ordering_genes = ordering_genes)
monocle_cds <-reduceDimension(monocle_cds,reduction_method = "DDRTree",max_components = 2)

monocle_cds <-orderCells(monocle_cds )
plot_cell_trajectory(monocle_cds, color_by = "State",cell_size = 0.75) 
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",cell_size = 0.75)
plot_cell_trajectory(monocle_cds, color_by = "group",cell_size = 0.75)
# remove(data,diff_test_res,fd,fData,pbmc.marker,pd,new.cluster.ids)
remove(pd,monocle_cds,fData,fd,data)
remove(ordering_genes,pbmcmarkers,pbmcmarkers_new,pbmcmarkers_new_new)
blast_genes <- row.names(subset(fData(monocle_cds),
                                gene_short_name %in% c("JMJD1C")))
plot_genes_jitter(monocle_cds[blast_genes,],
                  grouping = "State",
                  min_expr = 0.1)

my_genes <- row.names(subset(fData(monocle_cds),
                             gene_short_name %in% c("JMJD1C")))
cds_subset <- monocle_cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "group")
Tregtumor<-subset(pbmc_raw,idents = c("CD4_C12-CTLA4"))
Tregpbmc<-subset(pbmc_raw,idents = c("CD4_C10-FOXP3"))
Tregpbmc.data<-Tregpbmc@assays$RNA@counts
remove(Tregpbmc)
JMJD1C<-FetchData(Tregtumor,vars = "JMJD1C")
JMJD1C$sample<-rep("jmjd1c",length(rownames(JMJD1C)))

Namejmjd1chi<-rownames(JMJD1C[which(JMJD1C$JMJD1C>mean(JMJD1C$JMJD1C)),])
Namejmjd1clo<-rownames(JMJD1C[which(JMJD1C$JMJD1C<=mean(JMJD1C$JMJD1C)),])
JMJD1Chi.data<-as.data.frame(Tregtumor@assays$RNA@counts)[,Namejmjd1chi]
JMJD1Clo.data<-as.data.frame(Tregtumor@assays$RNA@counts)[,Namejmjd1clo]
remove(Namejmjd1chi,Namejmjd1clo,Tregtumor,JMJD1C)

JMJD1Chi.data <- as.data.frame(JMJD1Chi.data)
JMJD1Clo.data <- as.data.frame(JMJD1Clo.data)
Tregpbmc.data <- as.data.frame(Tregpbmc.data)

for (i in 1:length(colnames(JMJD1Chi.data))) {
  colnames(JMJD1Chi.data)[i] <- paste(colnames(JMJD1Chi.data)[i],"JMJD1Chi",i,sep = "-")  
}

for (i in 1:length(colnames(JMJD1Clo.data))) {
  colnames(JMJD1Clo.data)[i] <- paste(colnames(JMJD1Clo.data)[i],"JMJD1Clo",i,sep = "-")  
}

for (i in 1:length(colnames(Tregpbmc.data))) {
  colnames(Tregpbmc.data)[i] <- paste(colnames(Tregpbmc.data)[i],"Tregpbmc",i,sep = "-")  
}

JMJD1Chi.metadata<-data.frame(colnames(JMJD1Chi.data),rep("hi",length(colnames(JMJD1Chi.data))))
JMJD1Clo.metadata<-data.frame(colnames(JMJD1Clo.data),rep("lo",length(colnames(JMJD1Clo.data))))
Tregpbmc.metadata<-data.frame(colnames(Tregpbmc.data),rep("pbmc",length(colnames(Tregpbmc.data))))
colnames(JMJD1Chi.metadata)<-c("barcode","group")
colnames(JMJD1Clo.metadata)<-c("barcode","group")
colnames(Tregpbmc.metadata)<-c("barcode","group")

rownames(JMJD1Chi.metadata)<-JMJD1Chi.metadata[,1]
rownames(JMJD1Clo.metadata)<-JMJD1Clo.metadata[,1]
rownames(Tregpbmc.metadata)<-Tregpbmc.metadata[,1]


JMJD1Chi <- CreateSeuratObject(counts = JMJD1Chi.data, project = "IMMUNE_JMJD1Chi",meta.data = JMJD1Chi.metadata)
JMJD1Chi$type <- "JMJD1Chi"
JMJD1Chi <- NormalizeData(JMJD1Chi, verbose = FALSE)
JMJD1Chi <- FindVariableFeatures(JMJD1Chi, selection.method = "vst", nfeatures = 2000)


JMJD1Clo <- CreateSeuratObject(counts = JMJD1Clo.data, project = "IMMUNE_JMJD1Clo",meta.data = JMJD1Clo.metadata)
JMJD1Clo$type <- "JMJD1Clo"
JMJD1Clo <- NormalizeData(JMJD1Clo, verbose = FALSE)
JMJD1Clo <- FindVariableFeatures(JMJD1Clo, selection.method = "vst", nfeatures = 2000)

Tregpbmc <- CreateSeuratObject(counts = Tregpbmc.data, project = "IMMUNE_Tregpbmc",meta.data = Tregpbmc.metadata)
Tregpbmc$type <- "Tregpbmc"
Tregpbmc <- NormalizeData(Tregpbmc, verbose = FALSE)
Tregpbmc <- FindVariableFeatures(Tregpbmc, selection.method = "vst", nfeatures = 2000)
remove(JMJD1Chi.data,JMJD1Clo.data,JMJD1Chi.metadata,JMJD1Clo.metadata,Tregpbmc.data,Tregpbmc.metadata)
remove(i)

immune.anchors <- FindIntegrationAnchors(object.list = list(JMJD1Chi,JMJD1Clo,Tregpbmc), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- ScaleData(immune.combined, verbose = FALSE,assay = "RNA")
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.8)

p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "type")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
p1
p2
DimPlot(immune.combined, reduction = "umap", split.by = "type")

immune.combined <- RunTSNE(immune.combined, dims = 1:20)
DimPlot(immune.combined, reduction = "tsne",split.by = "type")
DimPlot(immune.combined, reduction = "tsne",group.by  = "type")
DimPlot(immune.combined, reduction = "tsne")

immune.combined<-RunLDA(immune.combined,labels = immune.combined@meta.data$group,assay = "RNA",features = rownames(immune.combined))
immune.combined<-RunLDA(immune.combined,labels = immune.combined@meta.data$group,assay = "integrated",features = rownames(immune.combined),reduction.name = "LDA_integrated")
immune.combined<-RunTSNE(immune.combined,reduction = "lda",reduction.name = "lda_tsne",dims = 1:2)
immune.combined<-RunTSNE(immune.combined,reduction = "LDA_integrated",reduction.name = "lda_tsne_integrated",dims = 1:2)
Idents(immune.combined)<-immune.combined$group
DimPlot(immune.combined,reduction = "lda_tsne")
DimPlot(immune.combined,reduction = "lda_tsne_integrated")
remove(Tregpbmc,p1,p2,JMJD1Chi,JMJD1Clo,immune.anchors)
Treghilo<-subset(pbmc,idents = c("hi","lo"))
write.table(Treghilo@assays$RNA@scale.data,"c:/Users/xjmik/Desktop/hilo.txt",sep = "\t")

Tregtumor<-subset(pbmc_raw,idents = c("CD4_C12-CTLA4"))
CD8<-subset(pbmc_raw,idents = c("CD8_C06-CD160","CD8_C04-GZMK","CD8_C05-CD6","CD8_C01-LEF1","CD8_C03-CX3CR1","CD8_C02-GPR183","CD8_C08-SLC4A10","CD8_C07-LAYN"))
Idents(CD8)<-CD8@meta.data$Bulksample
CD8tumor<-subset(CD8,idents = "T")
remove(CD8)
CD4<-subset(pbmc_raw,idents = c("CD4_C11-IL10","CD4_C05-CXCR6","CD4_C06-CXCR5","CD4_C04-TCF7","CD4_C02-ANXA1","CD4_C07-GZMK","CD4_C01-CCR7","CD4_C12-CTLA4","CD4_C03-GNLY","CD4_C10-FOXP3","CD4_C09-CXCL13","CD4_C08-IL23R"))
Idents(CD4)<-CD4@meta.data$Bulksample
CD4tumor<-subset(CD4,idents = "T")
remove(CD4)
Idents(Tregtumor)<-Tregtumor$Patient_ID
Idents(CD8tumor)<-CD8tumor$Patient_ID
Idents(CD4tumor)<-CD4tumor$Patient_ID
patient<-unique(Tregtumor@meta.data$Patient_ID)

for (i in 1:length(patient)) {
  assign(paste("Treg_",i,"_TIL",sep = ""),subset(Tregtumor,idents = patient[i]))
  assign("a",slot(get(paste("Treg_",i,"_TIL",sep = "")),"assays"))
  assign("b",as.data.frame(t(as.data.frame(a$RNA@data))))
  assign("c",select(b,one_of("JMJD1C")))
  assign(paste("pbmc_CD8T_",i,"_TIL",sep = ""),subset(CD8tumor,idents = patient[i]))
  assign("d",length(colnames(get(paste("pbmc_CD8T_",i,"_TIL",sep = ""))))/length(colnames(get(paste("Treg_",i,"_TIL",sep = "")))))
  assign(paste("scaldata_",i,"_TIL",sep = ""),data.frame(mean(c$JMJD1C),d))
  
}
scaldata<-scaldata_1_TIL
for (i in 2:length(patient)) {
  scaldata<-rbind(scaldata,get(paste("scaldata_",i,"_TIL",sep = "")))
}
colnames(scaldata)<-c("JMJD1C","CD8_Treg")
library(ggplot2)
library(ggpubr)
ggplot(data = scaldata,aes(x=JMJD1C,y=CD8_Treg)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=TRUE,size=1.5,color="red")+stat_cor(data = scaldata,method = "spearman") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
