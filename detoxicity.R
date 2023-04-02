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
d<-d[2:559]
GTEX<-read.table("c:/Users/xjmik/Downloads/rna_tissue_gtex.tsv",sep = "\t",header = TRUE)
GTEX_new<-GTEX[,2:4]
Gene<-intersect(d,unique(GTEX_new$Gene.name))
remove(d,GTEX)
a<-filter(GTEX_new,Gene.name %in% Gene[1])
b<-as.data.frame(t(data.frame(a$TPM)))
rownames(b)<-Gene[1]
colnames(b)<-a$Tissue
remove(a)
for ( i in 2:length(Gene)) {
 a<-filter(GTEX_new,Gene.name %in% Gene[i])
 c<-as.data.frame(t(data.frame(a$TPM)))
 rownames(c)<-Gene[i]
 colnames(c)<-a$Tissue
 remove(a)
 d<-intersect(colnames(b),colnames(c))
 b<-b[,d]
 c<-c[,d]
 b<-rbind(b,c)
 remove(c)
}

