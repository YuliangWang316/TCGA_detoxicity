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

library(cBioPortalData)
library(AnVIL)
library(GSVA)
library(TCGAbiolinks)

cbio <- cBioPortal()
studies <- getStudies(cbio, buildReport = TRUE)
lihc <- cBioDataPack("lihc_tcga", ask = FALSE)
lihc_genexpression<-as.data.frame(lihc@ExperimentList@listData$mrna_seq_v2_rsem_zscores_ref_diploid_samples@assays@data@listData,check.names =FALSE)
lihc_genexpression<-na.omit(lihc_genexpression)
lihc_clinical<-as.data.frame(lihc@colData@listData)
samplesTP <- TCGAquery_SampleTypes(colnames(lihc_genexpression), typesample = c("TP"))
lihc_genexpression<-lihc_genexpression[,samplesTP]
d<-intersect(d,rownames(lihc_genexpression))
lihc_genexpression_new<-lihc_genexpression[d,]

rownames(lihc_clinical)<-lihc_clinical$SAMPLE_ID
lihc_clinical<-lihc_clinical[samplesTP,]

lihc_genexpression_new<-as.data.frame(t(lihc_genexpression_new))
remove(lihc_genexpression)
sample<-intersect(rownames(lihc_clinical),rownames(lihc_genexpression_new))
lihc_clinical<-lihc_clinical[sample,]
lihc_genexpression_new<-lihc_genexpression_new[sample,]
lihc_genexpression_new<-na.omit(lihc_genexpression_new)
lihc_clinical<-cbind(lihc_clinical,lihc_genexpression_new)
remove(sample,lihc_genexpression_new)
remove(samplesTP,d)
setwd("c:/Users/xjmik/Desktop/Detoxicity_2/LIHC_tcga_OS/")
for (i in c(101:115,117:245,248:397,399:424,426:length(colnames(lihc_clinical)))) {
  df<-subset(lihc_clinical,select =c("OS_STATUS","OS_MONTHS","DFS_STATUS","DFS_MONTHS",colnames(lihc_clinical)[i]))
  
  df <- df[!is.na(df[,5]),]
  
  df<-df[which(df$OS_MONTHS!="[Not Available]"),]
  df<-df[which(df$OS_MONTHS!= 0),]
  
  
  for (j in 1:length(rownames(df))) {
    if(df$OS_STATUS[j] == "0:LIVING"){
      df$events[j]<-0
    }else if(df$OS_STATUS[j] == "1:DECEASED"){
      df$events[j]<-1
    }
  }
  remove(j)
  library(survival)
  library(survminer)
  df$OS_MONTHS<-as.numeric(df$OS_MONTHS)
  res.cut<-surv_cutpoint(df,time = "OS_MONTHS",event = "events",variables = colnames(df)[5])
  summary(res.cut)
  res.cat<-surv_categorize(res.cut)
  for (j in 1:length(rownames(res.cat))) {
    if(res.cat[,3][j] == "high"){
      res.cat[,3][j]<-TRUE
    }else if(res.cat[,3][j] == "low"){
      res.cat[,3][j]<-FALSE
    }
  }
  remove(j)
  cutoff<-res.cat[,3]
  fit<-survfit(Surv(OS_MONTHS,events)~ cutoff ,data = res.cat)
 
  p<-ggsurvplot(fit,data = res.cat,risk.table = TRUE,pval = TRUE,)
  pdf(file=paste(colnames(lihc_clinical)[i],"_OS.pdf",sep=""))
  print(p,newpage = FALSE)
  dev.off()
}
remove(df,fit,p,res.cat,res.cut,i,cutoff)
setwd("c:/Users/xjmik/Desktop/Detoxicity_2/LIHC_tcga_DFS/")

for (i in c(101:115,117:245,248:397,399:424,426:length(colnames(lihc_clinical)))) {
  df<-subset(lihc_clinical,select =c("OS_STATUS","OS_MONTHS","DFS_STATUS","DFS_MONTHS",colnames(lihc_clinical)[i]))
  
  df <- df[!is.na(df[,5]),]

  df<-df[which(df$DFS_MONTHS!="[Not Available]"),]
  df<-df[which(df$DFS_MONTHS!= 0),]


for (j in 1:length(rownames(df))) {
  if(df$DFS_STATUS[j] == "0:DiseaseFree"){
    df$events[j]<-0
  }else if(df$DFS_STATUS[j] == "1:Recurred/Progressed"){
    df$events[j]<-1
  }
}
  remove(j)
library(survival)
library(survminer)
df$DFS_MONTHS<-as.numeric(df$DFS_MONTHS)
res.cut<-surv_cutpoint(df,time = "DFS_MONTHS",event = "events",variables = colnames(df)[5])
summary(res.cut)
res.cat<-surv_categorize(res.cut)
for (j in 1:length(rownames(res.cat))) {
  if(res.cat[,3][j] == "high"){
    res.cat[,3][j]<-TRUE
  }else if(res.cat[,3][j] == "low"){
    res.cat[,3][j]<-FALSE
  }
}
remove(j)
cutoff<-res.cat[,3]
fit<-survfit(Surv(DFS_MONTHS,events)~ cutoff ,data = res.cat)

p<-ggsurvplot(fit,data = res.cat,risk.table = TRUE,pval = TRUE,)
pdf(file=paste(colnames(lihc_clinical)[i],"_DFS.pdf",sep=""))
print(p,newpage = FALSE)
dev.off()
}
remove(df,fit,p,res.cat,res.cut,i,cutoff)
