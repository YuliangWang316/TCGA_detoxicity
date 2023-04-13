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
GTEX<-read.table("c:/Users/xjmik/Downloads/rna_tissue_gtex.tsv",sep = "\t",header = TRUE)
GTEX_new<-GTEX[,c(2,3,6)]
Gene<-intersect(d,unique(GTEX_new$Gene.name))
remove(d,GTEX)
a<-filter(GTEX_new,Gene.name %in% Gene[1])
b<-as.data.frame(t(data.frame(a$nTPM)))
rownames(b)<-Gene[1]
colnames(b)<-a$Tissue
remove(a)
for ( i in 2:length(Gene)) {
  a<-filter(GTEX_new,Gene.name %in% Gene[i])
  c<-as.data.frame(t(data.frame(a$nTPM)))
  rownames(c)<-Gene[i]
  colnames(c)<-a$Tissue
  remove(a)
  d<-intersect(colnames(b),colnames(c))
  b<-b[,d]
  c<-c[,d]
  b<-rbind(b,c)
  remove(c)
}
remove(i,Gene,d)
b<-b[rowMeans(b)>0,]
write.table(b,file = "c:/Users/xjmik/Desktop/HPA_detoxicity.txt",sep = "\t")

library(UCSCXenaTools)
data(XenaData)
Sample<-XenaGenerate(subset = XenaHostNames=="gdcHub")
LIHC<-XenaFilter(Sample,filterCohorts = "GDC TCGA Liver Cancer") 
XenaFilter(LIHC,filterDatasets  = "TCGA-LIHC.htseq_fpkm.tsv") %>%
  XenaQuery() %>%
  XenaDownload() %>%
  XenaPrepare()-> LIHC_expression
XenaFilter(LIHC,filterDatasets  = "TCGA-LIHC.survival.tsv") %>%
  XenaQuery() %>%
  XenaDownload() %>%
  XenaPrepare()-> LIHC_Clinical
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

LIHC_expression<-as.data.frame(LIHC_expression)
rownames(LIHC_expression)<-LIHC_expression[,1]
LIHC_expression<-LIHC_expression[,-1]

LIHC_expression_gl<-TCGAanalyze_Normalization(LIHC_expression,geneInfoHT,method =  'geneLength')
genename<-rownames(LIHC_expression_gl)
e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
e<-e[!duplicated(e$SYMBOL),]
LIHC_expression_gl<-LIHC_expression_gl[e$ENSEMBL,]
rownames(LIHC_expression_gl)<-e$SYMBOL
remove(e,genename)
LIHC_expression_gl<-as.data.frame(LIHC_expression_gl)

newname_gl<-filter_tcga_barcodes(colnames(LIHC_expression_gl),analyte_target = "RNA")
LIHC_expression_gl_new<-LIHC_expression_gl[,newname_gl]
LIHC_expression_gl_new_new<-LIHC_expression_gl_new[rowMeans(LIHC_expression_gl_new)>0,]
