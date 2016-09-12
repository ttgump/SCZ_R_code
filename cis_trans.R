rm(list=ls())
library(data.table)
library(plyr)
#library(plyr)
#library(stringr)
setwd("~/Documents/Research/RNASeq-SCZ/Analysis/stringtieOutGencode/cuffdiff/summary/distance_correlation/")
chrs = c("chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
          "chr19", "chr2", "chr20", "chr21", "chr22", "chr3", "chr4", "chr5", "chr6", "chr7",
          "chr8", "chr9", "chrX", "chrY")

for(c in chrs) {
  if(!exists("distance.t")) {
    distance.t = read.csv(paste(c, "_lnc_coding_pair.csv", sep=""), stringsAsFactors = F)
  }
  else {
    temp = read.csv(paste(c, "_lnc_coding_pair.csv", sep=""), stringsAsFactors = F)
    distance.t = rbind(distance.t, temp)
  }
}
# distance.t = as.data.table(distance.t)
# setkey(distance.t,a,b)

load("lnc.cor.Rdata")

cor.table = melt(lnc.cor.table)
cor.table$distance = -1
cor.table = as.data.table(cor.table)
setkey(cor.table, Var1, Var2)

for(i in 1:nrow(distance.t)) {
  # cor.table[Var1 == distance.t[i,1] & Var2 == distance.t[i,2], distance := distance.t[i,3]]
  # cor.table[list(distance.t[i,1], distance.t[i,2])]$distance = distance.t[i,3]
  idx = cor.table[list(distance.t[i,1], distance.t[i,2]), which=T]
  cor.table[idx, distance := distance.t[i,3]]
}

write.csv(cor.table, "cor.table.csv", row.names=F)

dat = data.frame(cor=c(cis, trans), type=c(rep("cis", length(cis)), rep("trans", length(trans))))
tiff("cis_trans_cor.tiff", width=5, height=6, res=600,  units="in", compression="lzw", type="cairo")
ggplot(dat, aes(x=cor,colour=type)) + geom_density() + xlim(-1,1)
dev.off()