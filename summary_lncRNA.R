rm(list=ls())
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(stringr)
library(pheatmap)
library(openxlsx)
setwd("~/Documents/Research/RNASeq-SCZ/Analysis/stringtieOutGencode/NovellncRNAs/")

gtf.f <- import.gff("gencode.v24.novel.lncRNA.gtf")
lncRNA.gtf = as.data.frame(import.gff("~/Documents/Research/RNASeq-SCZ/hg38_annotation/gencode.v24.long_noncoding_RNAs.gtf"))
gtf <- as.data.frame(gtf.f)
exons <- gtf[gtf$type == 'exon',]

exons$feature_length <- exons$end - exons$start + 1
exons$exon_number <- 1
transcript_len<- aggregate(list(transcript_length= exons$feature_length), by= list(transcript_id= exons$transcript_id), sum)
exon_count <- aggregate(list(exon_number= exons$exon_number), by= list(transcript_id= exons$transcript_id), sum)

gene_name = gtf[!duplicated(gtf$transcript_id), ]
gene_name = gene_name[, c("transcript_id", "gene_type")]
gene_name = gene_name[!is.na(gene_name$transcript_id),]
gene_name$transcript_id = as.character(gene_name$transcript_id)
lncRNA.id = unique(lncRNA.gtf$transcript_id)
gene_name$gene_type[is.na(gene_name$gene_type)] = "lncRNA"
gene_name$gene_type[grepl("MSTRG", gene_name$transcript_id)] = "novel lncRNA"
gene_name$gene_type[gene_name$gene_type=="protein_coding"] = "protein coding gene"
gene_name = gene_name[gene_name$gene_type %in% c("protein coding gene", "novel lncRNA"),]
# write.csv(gene_name, "transcript_name.csv", row.names=F)

transcript_len = merge(transcript_len, gene_name, by.x=1, by.y=1, sort=F, suffixes="")
exon_count = merge(exon_count, gene_name, by.x=1, by.y=1, sort=F, suffixes="")
exon_count$gene_type = factor(exon_count$gene_type, levels=c("lncRNA", "novel lncRNA", "protein coding gene"),
                                  ordered =T)

exon_length = exons[,c("transcript_id", "feature_length")]
exon_length = merge(exon_length, gene_name, sort=F, suffixes="")
# exon_length = exon_length[exon_length$gene_type %in% c("protein_coding", "novel lncRNA"),]
exon_length$gene_type = factor(exon_length$gene_type, levels=c("lncRNA", "novel lncRNA", "protein coding gene"),
                               ordered =T)
transcript_len$gene_type = factor(transcript_len$gene_type, levels=c("lncRNA", "novel lncRNA", "protein coding gene"),
                               ordered =T)

tiff("SCZ_exon_length.tiff", width=5, height=3.5, res=600,  units="in", compression="lzw", type="cairo")
ggplot(exon_length, aes(x=feature_length, colour=gene_type)) + geom_density(na.rm=T) +
  ggtitle("Exon Size distributions") + xlim(1, 3000) +
  xlab("Size(bp)") + ylab("Density")
dev.off()

tiff("SCZ_exon_length3.tiff", width=4, height=4, res=600,  units="in", compression="lzw", type="cairo")
ggplot(exon_length, aes(x=gene_type, y=feature_length, fill=gene_type)) + geom_boxplot(na.rm=T, outlier.size=0.5) +
  ggtitle("Exon Size") + theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  ylab("Size(bp)") + ylim(0, 2000) + scale_fill_discrete(drop=FALSE)
dev.off()

median(exon_length[exon_length$gene_type=="protein coding gene", "feature_length"])
median(exon_length[exon_length$gene_type=="lncRNA", "feature_length"])
median(exon_length[exon_length$gene_type=="novel lncRNA", "feature_length"])
wilcox.test(x=exon_length[exon_length$gene_type=="protein_coding", "feature_length"],
        y=exon_length[exon_length$gene_type=="lncRNA", "feature_length"], alternative="less")
wilcox.test(x=exon_length[exon_length$gene_type=="protein_coding", "feature_length"],
        y=exon_length[exon_length$gene_type=="novel lncRNA", "feature_length"], alternative="less")

transcript_len = transcript_len[transcript_len$gene_type %in% c("protein_coding", "lncRNA", "novel lncRNA"),]
tiff("SCZ_transcript_length.tiff", width=5, height=3.5, res=600,  units="in", compression="lzw", type="cairo")
ggplot(transcript_len, aes(x=transcript_length, colour=gene_type)) + geom_density(na.rm=T) +
  ggtitle("Transcript size distributions") + xlim(1, 10000) +
  xlab("Size(bp)") + ylab("Density")
dev.off()

tiff("SCZ_transcript_length3.tiff", width=4, height=4, res=600,  units="in", compression="lzw", type="cairo")
ggplot(transcript_len, aes(x=gene_type, y=transcript_length, fill=gene_type)) + geom_boxplot(na.rm=T, outlier.size=0.5) +
  ggtitle("Transcript Size") + theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  ylab("Size(bp)") + ylim(0, 10000) + geom_boxplot(outlier.size=0.5) + scale_fill_discrete(drop=FALSE)
dev.off()

median(transcript_len[transcript_len$gene_type=="protein_coding", "transcript_length"])
median(transcript_len[transcript_len$gene_type=="lncRNA", "transcript_length"])
median(transcript_len[transcript_len$gene_type=="novel lncRNA", "transcript_length"])
wilcox.test(x=transcript_len[transcript_len$gene_type=="protein_coding", "transcript_length"],
        y=transcript_len[transcript_len$gene_type=="lncRNA", "transcript_length"], alternative="greater")
wilcox.test(x=transcript_len[transcript_len$gene_type=="protein_coding", "transcript_length"],
        y=transcript_len[transcript_len$gene_type=="novel lncRNA", "transcript_length"], alternative="greater")
wilcox.test(x=transcript_len[transcript_len$gene_type=="novel lncRNA", "transcript_length"],
            y=transcript_len[transcript_len$gene_type=="lncRNA", "transcript_length"], alternative="greater")

exon_count = exon_count[exon_count$gene_type %in% c("protein_coding", "lncRNA", "novel lncRNA"),]

tiff("SCZ_exon_number3.tiff", width=5.5, height=4, res=600,  units="in", compression="lzw", type="cairo")
ggplot(exon_count, aes(x=exon_number, fill=gene_type, y=..density.., width=0.85)) + 
  geom_histogram(binwidth=1, alpha=.5, position="identity") + 
  ggtitle("Exons per transcript") + xlim(1, 40) +
  xlab("Number of exons") + ylab("Density") + scale_fill_discrete(drop=FALSE)
dev.off()

median(exon_count[exon_count$gene_type=="protein_coding", "exon_number"])
median(exon_count[exon_count$gene_type=="lncRNA", "exon_number"])
median(exon_count[exon_count$gene_type=="novel lncRNA", "exon_number"])
wilcox.test(x=exon_count[exon_count$gene_type=="protein_coding", "exon_number"],
        y=exon_count[exon_count$gene_type=="lncRNA", "exon_number"])
wilcox.test(x=exon_count[exon_count$gene_type=="protein_coding", "exon_number"],
        y=exon_count[exon_count$gene_type=="novel lncRNA", "exon_number"])

orf.d = read.table("CPAT/gencode.v24.novel.lncRNA.txt", header=T, sep="\t", quote="")
orf.d = merge(orf.d, gene_name, by.x=0, by.y=1, sort=F, suffixes="")

orf.d = orf.d[orf.d$gene_type %in% c("protein coding gene", "lncRNA", "novel lncRNA"),]
orf.d$gene_type = factor(orf.d$gene_type, levels=c("lncRNA", "novel lncRNA", "protein coding gene"),
                        ordered =T)
tiff("SCZ_transcript_ORF.tiff", width=5, height=3.5, res=600,  units="in", compression="lzw", type="cairo")
ggplot(orf.d, aes(x=ORF_size, colour=gene_type)) + geom_density(na.rm=T) + 
  ggtitle("ORF size of transcripts") + xlim(0, 3000) +
  xlab("ORF size (bp)") + ylab("Density")
dev.off()

tiff("SCZ_transcript_ORF3.tiff", width=4, height=4, res=600,  units="in", compression="lzw", type="cairo")
ggplot(orf.d, aes(x=gene_type, y=ORF_size, fill=gene_type)) + geom_boxplot(na.rm=T, outlier.size=0.5) +
  ggtitle("ORF size of transcripts") + theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  ylab("Size(bp)") + ylim(0, 3000) + geom_boxplot(outlier.size=0.5) + scale_fill_discrete(drop=FALSE)
dev.off()

median(orf.d[orf.d$gene_type=="protein_coding", "ORF_size"])
median(orf.d[orf.d$gene_type=="lncRNA", "ORF_size"])
median(orf.d[orf.d$gene_type=="novel lncRNA", "ORF_size"])
wilcox.test(x=orf.d[orf.d$gene_type=="protein_coding", "ORF_size"],
        y=orf.d[orf.d$gene_type=="lncRNA", "ORF_size"], alternative="greater")
wilcox.test(x=orf.d[orf.d$gene_type=="protein_coding", "ORF_size"],
        y=orf.d[orf.d$gene_type=="novel lncRNA", "ORF_size"], alternative="greater")

tiff("SCZ_transcript_Coding_Score3.tiff", width=4, height=4, res=600,  units="in", compression="lzw", type="cairo")
ggplot(orf.d, aes(x=gene_type, y=coding_prob, fill=gene_type)) + geom_boxplot(outlier.size=0.5) + 
  ggtitle("Coding potential of transcripts") + theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  ylab("CPAT score")
dev.off()

median(orf.d[orf.d$gene_type=="protein_coding", "coding_prob"])
median(orf.d[orf.d$gene_type=="lncRNA", "coding_prob"])
median(orf.d[orf.d$gene_type=="novel lncRNA", "coding_prob"])
wilcox.test(x=orf.d[orf.d$gene_type=="protein_coding", "coding_prob"],
            y=orf.d[orf.d$gene_type=="lncRNA", "coding_prob"], alternative="greater")
wilcox.test(x=orf.d[orf.d$gene_type=="protein_coding", "coding_prob"],
            y=orf.d[orf.d$gene_type=="novel lncRNA", "coding_prob"], alternative="greater")

fpkm = read.csv("../cuffdiff/summary/gene_fpkm.csv", row.names=1)
gene.name.gene = read.csv("../cuffdiff/gene_name.csv")
lncRNA.id.gene = read.csv("../../../hg38_annotation/lncRNA_gene_id.csv")$x
fpkm.name = merge(fpkm, gene.name.gene, by.x=0, by.y=1, sort=F)
fpkm.name$gene_type = as.character(fpkm.name$gene_type)
fpkm.name$gene_type[fpkm.name$Row.names %in% lncRNA.id.gene] = "lncRNA"
fpkm.name$gene_type[grepl("MSTRG", fpkm.name$Row.names)] = "novel lncRNA"
fpkm.name = fpkm.name[fpkm.name$gene_type %in% c("protein_coding", "lncRNA", "novel lncRNA"),]
fpkm.name$gene_type[fpkm.name$gene_type=="protein_coding"] = "protein coding gene"

gene_exp = data.frame(gene=rep(fpkm.name$Row.names, ncol(fpkm.name[,2:47])),
                      type=rep(fpkm.name$gene_type, ncol(fpkm.name[,2:47])),
                      FPKM=c(apply(fpkm.name[,2:47],2,rbind)),
                      log10.FPKM=log10(c(apply(fpkm.name[,2:47],2,rbind))+1e-5))
gene_exp1 = gene_exp[gene_exp$FPKM>0,]
gene_exp1$type = factor(gene_exp1$type, levels=c("lncRNA", "novel lncRNA", "protein coding gene"),
                        ordered =T)

median(gene_exp1[gene_exp1$type=="protein coding gene", "FPKM"])
median(gene_exp1[gene_exp1$type=="lncRNA", "FPKM"])
median(gene_exp1[gene_exp1$type=="novel lncRNA", "FPKM"])
wilcox.test(x=gene_exp1[gene_exp1$type=="protein_coding", "FPKM"],
            y=gene_exp1[gene_exp1$type=="lncRNA", "FPKM"], alternative="greater")
wilcox.test(x=gene_exp1[gene_exp1$type=="protein_coding", "FPKM"],
            y=gene_exp1[gene_exp1$type=="novel lncRNA", "FPKM"], alternative="greater")

tiff("SCZ_gene_FPKM.tiff", width=5, height=3.5, res=600,  units="in", compression="lzw", type="cairo")
ggplot(gene_exp1, aes(x=log10.FPKM, colour=type)) + geom_density(na.rm=T) + 
  ggtitle("Gene Expression distributions") + xlim(-4,4) +
  xlab("Expression (in log10 FPKM)") + ylab("Density")
dev.off()

tiff("SCZ_gene_FPKM3.tiff", width=4, height=4, res=600,  units="in", compression="lzw", type="cairo")
ggplot(gene_exp1, aes(x=type, y=log10.FPKM, fill=type)) + geom_boxplot(na.rm=T, outlier.size=0.5) +
  ggtitle("Gene Expression Levels") + theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  ylab("Expression (in log10 FPKM)") + ylim(-4,4) + geom_boxplot(outlier.size=0.5) +
  scale_fill_discrete(drop=FALSE)
dev.off()
