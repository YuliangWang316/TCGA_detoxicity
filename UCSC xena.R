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
  genset<-read.table("c:/Users/xjmik/Downloads/1C_signature_182.txt",sep = "\t")
colnames(genset)<-"J1C_signature"
geneset_list<-as.list(genset)

gsva_gl<-gsva(expr = as.matrix(LIHC_expression_gl_new_new),gset.idx.list = geneset_list,kcdf="Poisson",parallel.sz=20)
gsva_gl<-as.data.frame(t(gsva_gl))
gsva_gl$sample <- sapply(strsplit(rownames(gsva_gl),'-'),function(x) paste0(x[1:3],collapse="-"))
samplesTP <- TCGAquery_SampleTypes(colnames(LIHC_expression_gl_new_new), typesample = c("TP"))
rownames(LIHC_Clinical)<-LIHC_Clinical$sample
LIHC_Clinical<-LIHC_Clinical[samplesTP,]
colnames(LIHC_Clinical)[3]<-"PATIENT"
LIHC_Clinical$"J1C_signature" <- gsva_gl[match(LIHC_Clinical$PATIENT,gsva_gl$sample),][,"J1C_signature"]


df<-subset(LIHC_Clinical,select =c("OS","OS.time","J1C_signature"))

df <- df[!is.na(df$J1C_signature),]


library(survival)
library(survminer)

res.cut<-surv_cutpoint(df,time = "OS.time",event = "OS",variables = c("J1C_signature"))
summary(res.cut)
res.cat<-surv_categorize(res.cut)
fit<-survfit(Surv(OS.time,OS)~ J1C_signature ,data = res.cat)
ggsurvplot(fit,data = res.cat,risk.table = TRUE,pval = TRUE)
