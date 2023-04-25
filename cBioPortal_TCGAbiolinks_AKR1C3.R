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
lihc <- cBioDataPack("lihc_tcga_pan_can_atlas_2018", ask = FALSE)
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

sample<-intersect(rownames(lihc_clinical),rownames(lihc_genexpression_new))
lihc_clinical<-lihc_clinical[sample,]
lihc_genexpression_new<-lihc_genexpression_new[sample,]
lihc_genexpression_new<-na.omit(lihc_genexpression_new)
lihc_clinical<-cbind(lihc_clinical,lihc_genexpression_new)
remove(sample,lihc_genexpression_new)
remove(samplesTP,d)

df<-subset(lihc_clinical,select =c("OS_STATUS","OS_MONTHS","DFS_STATUS","DFS_MONTHS","AKR1C3","PFS_STATUS","PFS_MONTHS","DSS_STATUS","DSS_MONTHS"))

df <- df[!is.na(df$AKR1C3),]

df<-df[which(df$OS_MONTHS!="[Not Available]"),]
df<-df[which(df$OS_MONTHS!= 0),]


for (j in 1:length(rownames(df))) {
  if(df$OS_STATUS[j] == "0:LIVING"){
    df$events[j]<-0
  }else if(df$OS_STATUS[j] == "1:DECEASED"){
    df$events[j]<-1
  }
}
library(survival)
library(survminer)
df$OS_MONTHS<-as.numeric(df$OS_MONTHS)
res.cut<-surv_cutpoint(df,time = "OS_MONTHS",event = "events",variables = c("AKR1C3"))
summary(res.cut)
res.cat<-surv_categorize(res.cut)
fit<-survfit(Surv(OS_MONTHS,events)~ AKR1C3 ,data = res.cat)
ggsurvplot(fit,data = res.cat,risk.table = TRUE,pval = TRUE)

remove(cbio,j,fit,lihc,res.cut,studies,df)

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
hi<-rownames(res.cat[which(res.cat$AKR1C3 == "high"),])
lo<-rownames(res.cat[which(res.cat$AKR1C3 == "low"),])
remove(samplesTP)
remove(res.cat,count_gl_new_new)
b<-intersect(rownames(count_gl_new_new_tumor),rownames(lihc_genexpression))
count_gl_new_new_tumor_new<-count_gl_new_new_tumor[b,]
count_gl_new_new_tumor_new_hi<-as.data.frame(t(as.data.frame(t(count_gl_new_new_tumor_new))[hi,]))
count_gl_new_new_tumor_new_low<-as.data.frame(t(as.data.frame(t(count_gl_new_new_tumor_new))[lo,]))
count<-cbind(count_gl_new_new_tumor_new_hi,count_gl_new_new_tumor_new_low)
remove(count_gl_new_new_tumor_new,count_gl_new_new_tumor,count_gl_new_new_tumor_new_hi,count_gl_new_new_tumor_new_low)
lihc_genexpression_new<-lihc_genexpression[b,]
remove(b,lihc_genexpression)
lihc_genexpression_new_hi<-as.data.frame(t(as.data.frame(t(lihc_genexpression_new))[hi,]))
lihc_genexpression_new_lo<-as.data.frame(t(as.data.frame(t(lihc_genexpression_new))[lo,]))
lihc_gene<-cbind(lihc_genexpression_new_hi,lihc_genexpression_new_lo)
remove(hi,lo,lihc_genexpression_new,lihc_genexpression_new_hi,lihc_genexpression_new_lo)
condition<-factor(c(rep("hi",50),rep("lo",310)),levels = c("lo","hi"))
library(DESeq2)
colData<-data.frame(row.names = colnames(count),condition)

dds <- DESeqDataSetFromMatrix(count, colData, design= ~ condition)
dds <- DESeq(dds)

res= results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)
write.csv(res,file="c:/Users/xjmik/Desktop/Detoxicity_2/AKR1C3_survival_hilo.csv")
remove(colData,dds,res,condition)
write.table(lihc_gene,"c:/Users/xjmik/Desktop/Detoxicity_2/AKR1C3_survival_hilo_forGSEA.txt",sep = "\t")
v <-voom(count, plot = F, save.plot = F)
out=v$E
out=rbind(ID=colnames(out),out)
write.table(out,file="AKR1C3survivalhilo_count.txt",sep="\t",quote=F,col.names=F)  
source("TMEimmune31.CIBERSORT.R")
results=CIBERSORT("ref.txt", "AKR1C3survivalhilo_count.txt", perm=100, QN=TRUE)
remove(v,out,results)
results=CIBERSORT("ref.txt", "AKR1C3survivalhilo_cBioPortalgene.txt", perm=100, QN=TRUE)
remove(v,out,results)
