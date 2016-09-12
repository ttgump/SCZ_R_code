rm(list=ls())
library(GenomicRanges)
library(rtracklayer)
library(plyr)

setwd("~/Documents/Research/RNASeq-SCZ/Analysis/stringtieOutGencode/cuffdiff/summary/distance_correlation/")
gtf.t.0.0 = read.csv("gene_coordinate.csv", stringsAsFactors=F)
rownames(gtf.t.0.0) = gtf.t.0.0$gene_id
distGene <- function(x="", y="", data=NULL) {
  x.t = data[x,c("seqnames", "start","end")]
  y.t = data[y,c("seqnames", "start","end")]
  gr1 = GRanges(seqnames=x.t$seqnames, ranges=IRanges(start=x.t$start, end=x.t$end, names=x.t$gene_id))
  gr2 = GRanges(seqnames=y.t$seqnames, ranges=IRanges(start=y.t$start, end=y.t$end, names=y.t$gene_id))
  return(distance(gr1, gr2))
}

# load("coding.cor.Rdata")
load("lnc.cor.Rdata")

for(chr in sort(unique(gtf.t.0.0$seqnames))) {
  fileConn<-file(paste(chr, "R", sep="."))
  line =
paste("
rm(list=ls())
library(GenomicRanges)
library(rtracklayer)
library(plyr)
gtf.t.0.0 = read.csv(\"gene_coordinate.csv\", stringsAsFactors=F)
rownames(gtf.t.0.0) = gtf.t.0.0$gene_id
  distGene <- function(x=\"\", y=\"\", data=NULL) {
  x.t = data[x,c(\"seqnames\", \"start\",\"end\")]
  y.t = data[y,c(\"seqnames\", \"start\",\"end\")]
  gr1 = GRanges(seqnames=x.t$seqnames, ranges=IRanges(start=x.t$start, end=x.t$end, names=x.t$gene_id))
  gr2 = GRanges(seqnames=y.t$seqnames, ranges=IRanges(start=y.t$start, end=y.t$end, names=y.t$gene_id))
  return(distance(gr1, gr2))
  }
load(\"lnc.cor.Rdata\")
lncRNA.id1 = rownames(lnc.cor.table)[rownames(lnc.cor.table) %in% gtf.t.0.0[gtf.t.0.0$seqnames==\"",chr,"\", \"gene_id\"]]
coding.id1 = colnames(lnc.cor.table)[colnames(lnc.cor.table) %in% gtf.t.0.0[gtf.t.0.0$seqnames==\"",chr,"\", \"gene_id\"]]
lnc.coding.pair = arrange(expand.grid(a=lncRNA.id1, b=coding.id1, stringsAsFactors=F), a)
lnc.coding.pair$distance <- apply(lnc.coding.pair, 1 ,function(z) {
  distGene(x=as.character(z[1]), y=as.character(z[2]), data=gtf.t.0.0)
})
write.csv(lnc.coding.pair, paste(\"",chr,"\", \"lnc_coding_pair.csv\", sep=\"_\"), row.names=F)", sep="")
  writeLines(line, fileConn)
  close(fileConn)
}
# lncRNA.id1 = rownames(lnc.cor.table)
# coding.id1 = colnames(lnc.cor.table)
# lnc.coding.pair = expand.grid(lncRNA.id1, coding.id1)
# lnc.coding.pair$distance <- apply(lnc.coding.pair, 1 ,function(z) {
#   distGene(x=as.character(z[1]), y=as.character(z[2]), data=gtf.t.0.0)
# })
# lnc.coding.pair = lnc.coding.pair[complete(lnc.coding.pair),]
# 
# write.csv(lnc.coding.pair, "lnc_coding_pair.csv", row.names=F)