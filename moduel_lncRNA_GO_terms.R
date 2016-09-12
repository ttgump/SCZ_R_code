rm(list=ls())
library(stringr)
library(pheatmap)
library(NbClust)
library(openxlsx)
library(ggplot2)
setwd("~/Documents/Research/RNASeq-SCZ/Analysis/stringtieOutGencode/cuffdiff/lncRNA_GSEA")
enrich = read.csv("lncRNA_go_enrich_0.05.csv", stringsAsFactors=F)
enrich.h = enrich
rownames(enrich.h) = enrich$NAME
enrich.h = enrich.h[,-1]
# enrich.h = enrich.h[,colnames(enrich.h) %in% all.lncRNA.id]
enrich.h = enrich.h[rowSums(abs(enrich.h))>0, colSums(abs(enrich.h))>0]
enrich.h = as.data.frame(enrich.h)

ENSG00000225465.8 = enrich.h[,"ENSG00000225465.8",drop=F]
ENSG00000225465.8 = ENSG00000225465.8[ENSG00000225465.8!=0,,drop=F]

ENSG00000276269.1 = enrich.h[,"ENSG00000276269.1",drop=F]
ENSG00000276269.1 = ENSG00000276269.1[ENSG00000276269.1!=0,,drop=F]

ENSG00000280237.1 = enrich.h[,"ENSG00000280237.1",drop=F]
ENSG00000280237.1 = ENSG00000280237.1[ENSG00000280237.1!=0,,drop=F]

write.xlsx(list(ENSG00000225465.8=ENSG00000225465.8,
                ENSG00000276269.1=ENSG00000276269.1,
                ENSG00000280237.1=ENSG00000280237.1), "Module_lncRNA_enrich_GO_terms0.xlsx", row.names=T)

ENSG00000225465.8 = enrich.h[,"ENSG00000225465.8",drop=F]
ENSG00000225465.8 = ENSG00000225465.8[ENSG00000225465.8!=0,,drop=F]

ENSG00000267519.3 = enrich.h[,"ENSG00000267519.3",drop=F]
ENSG00000267519.3 = ENSG00000267519.3[ENSG00000267519.3!=0,,drop=F]

ENSG00000280237.1 = enrich.h[,"ENSG00000280237.1",drop=F]
ENSG00000280237.1 = ENSG00000280237.1[ENSG00000280237.1!=0,,drop=F]
ENSG00000225953.2 = enrich.h[,"ENSG00000225953.2",drop=F]
ENSG00000225953.2 = ENSG00000225953.2[ENSG00000225953.2!=0,,drop=F]
ENSG00000273409.1 = enrich.h[,"ENSG00000273409.1",drop=F]
ENSG00000273409.1 = ENSG00000273409.1[ENSG00000273409.1!=0,,drop=F]
ENSG00000228058.1 = enrich.h[,"ENSG00000228058.1",drop=F]
ENSG00000228058.1 = ENSG00000228058.1[ENSG00000228058.1!=0,,drop=F]
MSTRG.33355.1 = enrich.h[,"MSTRG.33355.1",drop=F]
MSTRG.33355.1 = MSTRG.33355.1[MSTRG.33355.1!=0,,drop=F]

write.xlsx(list(ENSG00000225465.8=ENSG00000225465.8,
                     ENSG00000267519.3=ENSG00000267519.3,
                     ENSG00000280237.1=ENSG00000280237.1,
                     ENSG00000225953.2=ENSG00000225953.2,
                     ENSG00000273409.1=ENSG00000273409.1,
                     ENSG00000228058.1=ENSG00000228058.1,
                     MSTRG.33355.1=MSTRG.33355.1), "Module_lncRNA_enrich_GO_terms.xlsx", row.names=T)