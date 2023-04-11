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
i<-"LIHC"
cancer_type=paste("TCGA",i,sep="-")

query <- GDCquery(project = cancer_type, 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")

GDCdownload(query, method = "api", files.per.chunk = 100)
expdat <- GDCprepare(query = query)
count_matrix<-as.data.frame(assay(expdat))
count_gl<-TCGAanalyze_Normalization(count_matrix, geneInfoHT,method =  'geneLength')
remove(count_matrix,expdat,query)
genename<-rownames(count_gl)
e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
e<-e[!duplicated(e$SYMBOL),]
count_gl<-count_gl[e$ENSEMBL,]
rownames(count_gl)<-e$SYMBOL
remove(e,genename)
count_gl<-as.data.frame(count_gl)

newname_gl<-filter_tcga_barcodes(colnames(count_gl),analyte_target = "RNA")
count_gl_new<-count_gl[,newname_gl]
count_gl_new_new<-count_gl_new[rowMeans(count_gl_new)>0,]
remove(count_gl,count_gl_new,cancer_type,i,newname_gl,request_cancer)

samplesTP <- TCGAquery_SampleTypes(colnames(count_gl_new_new), typesample = c("TP"))
count_gl_new_new_tumor<-count_gl_new_new[,samplesTP]

samplesNT <- TCGAquery_SampleTypes(colnames(count_gl_new_new), typesample = c("NT"))
count_gl_new_new_normal<-count_gl_new_new[,samplesNT]
remove(samplesNT,samplesTP)
remove(count_gl_new_new)

genename_1<-rownames(count_gl_new_new_normal)

metadata_normal<-data.frame(colnames(count_gl_new_new_normal),rep("LIHC",length(colnames(count_gl_new_new_normal))))
rownames(metadata_normal)<-metadata_normal[,1]
colnames(metadata_normal)<-c("Barcode","group")

metadata_tumor<-data.frame(colnames(count_gl_new_new_tumor),rep("LIHC",length(colnames(count_gl_new_new_tumor))))
rownames(metadata_tumor)<-metadata_tumor[,1]
colnames(metadata_tumor)<-c("Barcode","group")


normal<- CreateSeuratObject(counts = count_gl_new_new_normal, project = "pbmc3k",meta.data = metadata_normal)
Tumor<- CreateSeuratObject(counts = count_gl_new_new_tumor, project = "pbmc3k", meta.data = metadata_tumor)  

normal <- NormalizeData(normal)  
Tumor <- NormalizeData(Tumor)  

normal <- FindVariableFeatures(normal, selection.method = "vst", nfeatures = 2000)
Tumor <- FindVariableFeatures(Tumor, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(normal)
normal <- ScaleData(normal, features = all.genes)  
Tumor <- ScaleData(Tumor, features = all.genes)  
remove(count_gl_new_new_normal,count_gl_new_new_tumor,metadata_normal,metadata_tumor,all.genes)  

normal <- RunPCA(normal, features = VariableFeatures(object = normal))  
Tumor <- RunPCA(Tumor, features = VariableFeatures(object = Tumor))  

Idents(normal)<-normal$group
Idents(Tumor)<-Tumor$group

LIHC_tumor<-subset(Tumor,idents = "LIHC")
LIHC_normal<-subset(normal,idents = "LIHC")

LIHC_count<-cbind(LIHC_tumor@assays$RNA@counts,LIHC_normal@assays$RNA@counts)  
LIHC_normal$GROUP<-rep("normal",length(colnames(LIHC_normal)))
LIHC_tumor$GROUP<-rep("tumor",length(colnames(LIHC_tumor)))

LIHC_metadata<-rbind(LIHC_tumor@meta.data[,c(1:6)],LIHC_normal@meta.data[,c(1:6)])  

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

gene<-intersect(d,rownames(pbmc))
remove(genename_1)
pbmc.markers<-FindAllMarkers(pbmc,only.pos = TRUE,min.pct = 0,logfc.threshold = 0)  
pbmc.tumor_vs_normal<-FindMarkers(pbmc,min.pct = 0,logfc.threshold = 0,ident.1 = "tumor",ident.2 = "normal")  
pbmc.markers_detoxicity<-FindAllMarkers(pbmc,only.pos = TRUE,min.pct = 0,logfc.threshold = 0,features = gene)  
pbmc.tumor_vs_normal_detoxicity<-FindMarkers(pbmc,min.pct = 0,logfc.threshold = 0,ident.1 = "tumor",ident.2 = "normal",features = gene)  
remove(all.genes)

write.table(pbmc.markers,file = "c:/Users/xjmik/Desktop/LIHC_Allmarker_tumor_normal.txt",sep = "\t")
write.table(pbmc.markers_detoxicity,file = "c:/Users/xjmik/Desktop/LIHC_Allmarker_tumor_normal_detoxicity.txt",sep = "\t")
write.table(pbmc.tumor_vs_normal,file = "c:/Users/xjmik/Desktop/LIHC_tumor_vs_normal.txt",sep = "\t")
write.table(pbmc.tumor_vs_normal_detoxicity,file = "c:/Users/xjmik/Desktop/LIHC_tumor_vs_normal_detoxicity.txt",sep = "\t")
remove(pbmc.markers,pbmc.markers_detoxicity,pbmc.tumor_vs_normal,pbmc.tumor_vs_normal_detoxicity)

gene_2<-c("AKR1C3","GPX4","HTATIP2","GPX1","FASN","AKR1B10","ACSL4","PON2")


setwd("D:/TEST3/")
v <-voom(Tumor@assays$RNA@counts, plot = F, save.plot = F)
out=v$E
out=rbind(ID=colnames(out),out)
write.table(out,file="uniq.symbol_LIHC_tumor.txt",sep="\t",quote=F,col.names=F)  
source("TMEimmune31.CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol_LIHC_tumor.txt", perm=100, QN=TRUE)
remove(v,out,results)

v <-voom(normal@assays$RNA@counts, plot = F, save.plot = F)
out=v$E
out=rbind(ID=colnames(out),out)
write.table(out,file="uniq.symbol_LIHC_normal.txt",sep="\t",quote=F,col.names=F)  
source("TMEimmune31.CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol_LIHC_normal.txt", perm=100, QN=TRUE)
remove(v,out,results)

v <-voom(cbind(as.data.frame(Tumor@assays$RNA@counts),as.data.frame(normal@assays$RNA@counts)), plot = F, save.plot = F)
out=v$E
out=rbind(ID=colnames(out),out)
write.table(out,file="uniq.symbol_LIHC_Tumor_normal.txt",sep="\t",quote=F,col.names=F)  
source("TMEimmune31.CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol_LIHC_Tumor_normal.txt", perm=100, QN=TRUE)
remove(v,out,results)
remove(d,pbmc)

gene3<-c("AKR1C3","HTATIP2","AKR1B10","ACSL4")
Tumor_Cibersort<-read.table("c:/Users/xjmik/Desktop/Detoxicity_2/LIHC_tumor_CIBERSORT-Results.txt",sep = "\t",header = TRUE,row.names = 1)
Tumor_Cibersort<-Tumor_Cibersort[1:22]
Tumor_genecount<-as.data.frame(Tumor@assays$RNA@data)[gene,]
Tumor_genecount<-Tumor_genecount[gene3,]
Tumor_genecount<-as.data.frame(t(Tumor_genecount))

library(psych)
result<-corr.test(Tumor_Cibersort,Tumor_genecount,use = "pairwise",method = "spearman",adjust = "fdr",alpha = .05,ci = TRUE,minlength = 10)
result<-as.data.frame(print(result,short = FALSE))
result<-result[which(result$raw.p < 0.05),]
write.table(result,file = "c:/Users/xjmik/Desktop/Detoxicity_2/Tumor_immmune_corr_gene.txt",sep = "\t")
remove(Tumor_Cibersort,Tumor_genecount,result)

Normal_Cibersort<-read.table("c:/Users/xjmik/Desktop/Detoxicity_2/LIHC_Normal_CIBERSORT-Results.txt",sep = "\t",header = TRUE,row.names = 1)
Normal_Cibersort<-Normal_Cibersort[1:22]
Normal_genecount<-as.data.frame(normal@assays$RNA@data)[gene,]
Normal_genecount<-Normal_genecount[gene3,]
Normal_genecount<-as.data.frame(t(Normal_genecount))

library(psych)
result<-corr.test(Normal_Cibersort,Normal_genecount,use = "pairwise",method = "spearman",adjust = "fdr",alpha = .05,ci = TRUE,minlength = 10)
result<-as.data.frame(print(result,short = FALSE))
result<-result[which(result$raw.p < 0.05),]
write.table(result,file = "c:/Users/xjmik/Desktop/Detoxicity_2/Normal_immmune_corr_gene.txt",sep = "\t")
remove(Normal_Cibersort,Normal_genecount,result)
