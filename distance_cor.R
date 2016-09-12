rm(list=ls())
library(GenomicRanges)
library(rtracklayer)
library(plyr)

setwd("~/Documents/Research/RNASeq-SCZ/Analysis/stringtieOutGencode/cuffdiff/summary/distance_correlation")
gtf <- import.gff("~/Documents/Research/RNASeq-SCZ/Analysis/stringtieOutGencode/NovellncRNAs/gencode.v24.novel.lncRNA.gtf", format="gtf")
gtf.t <- as.data.frame(gtf)
gtf.t1 <- gtf.t[gtf.t$type=="exon",]
gtf.t1 <- gtf.t1[!is.na(gtf.t1$gene_id),]
gtf.t0 <- gtf.t1[,c("gene_id", "seqnames", "start","end")]
gtf.t0.1 <- aggregate(gtf.t0$start, by=list(gtf.t0$gene_id), FUN=min)
colnames(gtf.t0.1) <- c("gene_id", "start")
gtf.t0.2 <- aggregate(gtf.t0$end, by=list(gtf.t0$gene_id), FUN=max)
colnames(gtf.t0.2) <- c("gene_id", "end")
gtf.t.0.0 <- merge(gtf.t0.1, gtf.t0.2, sort=F, surffix="")
gtf.tmp <- gtf.t0[,c("gene_id","seqnames")]
gtf.tmp <- gtf.tmp[!duplicated(gtf.tmp$gene_id),]
gtf.t.0.0 <- merge(gtf.t.0.0, gtf.tmp, by="gene_id", sort=F, all.x=T, all.y=F, surffix="")
write.csv(gtf.t.0.0, "gene_coordinate.csv", row.names=F)

gtf.t.0.0 = read.csv("gene_coordinate.csv")
rownames(gtf.t.0.0) = gtf.t.0.0$gene_id
distGene <- function(x="", y="", data=NULL) {
  x.t = data[x,c("seqnames", "start","end")]
  y.t = data[y,c("seqnames", "start","end")]
  gr1 = GRanges(seqnames=x.t$seqnames, ranges=IRanges(start=x.t$start, end=x.t$end, names=x.t$gene_id))
  gr2 = GRanges(seqnames=y.t$seqnames, ranges=IRanges(start=y.t$start, end=y.t$end, names=y.t$gene_id))
  return(distance(gr1, gr2))
}

gene.name = read.csv("../../gene_name.csv")

lncRNA.id = as.vector(read.csv("~/Documents/Research/RNASeq-SCZ/hg38_annotation/lncRNA_gene_id.csv")$x)
coding.id = as.vector(gene.name[gene.name$gene_type=="protein_coding", "gene_id"])


setwd("~/Documents/Research/RNASeq-SCZ/Analysis/stringtieOutGencode/cuffdiff/summary/")
fpkm = read.csv("gene_fpkm.csv", row.names=1)
# fpkm = fpkm[,-which(colnames(fpkm) %in% c("Sample_R10673"))]

All.expressed_gene = rownames(fpkm)[apply(fpkm, 1, function(z) {quantile(z, 0.95)>0.5 & quantile(z, 0.5)>0})]
All.expressed_coding = All.expressed_gene[All.expressed_gene %in% gene.name[gene.name$gene_type=="protein_coding", "gene_id"]]
All.expressed_lncRNAs = All.expressed_gene[All.expressed_gene %in% lncRNA.id]

coding.fpkm = fpkm[rownames(fpkm) %in% All.expressed_coding,]
coding.fpkm1 = merge(coding.fpkm, gene.name, by.x=0, by.y="gene_id")
coding.fpkm1$geneMean = rowMeans(coding.fpkm1[,-c(1,48:49)])
coding.fpkm1 = coding.fpkm1[order(coding.fpkm1$geneMean, decreasing=T),]
coding.fpkm1 = coding.fpkm1[!duplicated(coding.fpkm1$gene_name),]
rownames(coding.fpkm1) = coding.fpkm1$Row.names
coding.fpkm1 = coding.fpkm1[,-c(1,48:50)]
coding.fpkm1 = log2(coding.fpkm1+1)
lncRNA.fpkm = fpkm[rownames(fpkm) %in% All.expressed_lncRNAs,]
lncRNA.fpkm = log2(lncRNA.fpkm+1)
lnc.cor.table = cor(t(lncRNA.fpkm), t(coding.fpkm1), method="pearson")
coding.cor.table = cor(t(coding.fpkm1), method="pearson")

save(lnc.cor.table, file="lnc.cor.Rdata")
save(coding.cor.table, file="coding.cor.Rdata")

setwd("~/Documents/Research/RNASeq-SCZ/Analysis/stringtieOutGencode/cuffdiff/summary/distance_correlation")
lncRNA.id1 = rownames(lncRNA.fpkm)
coding.id1 = rownames(coding.fpkm)
lnc.coding.pair = expand.grid(lncRNA.id1, coding.id1)
lnc.coding.pair$distance <- apply(lnc.coding.pair, 1 ,function(z) {
  distGene(x=as.character(z[1]), y=as.character(z[2]), data=gtf.t.0.0)
})

coding.pair = expand.grid(coding.id1, coding.id1)
coding.pair$distance <- apply(coding.pair, 2, function(z) {
  distGene(x=as.character(z[1]), y=as.character(z[2]), data=gtf.t.0.0)
})