#!/usr/env/bin/Rscript
# _*_ coding: utf-8 _*_
library('getopt')
 
command=matrix(c( 
  'help', 'h', 0,'loical', '帮助文档',
  'pdfinput', 'i', 1, 'character', '判断值的结果',
  'noloutput', 'o', 1, 'character', '标准化的判断值的结果'),byrow=TRUE,ncol=5)
args=getopt(command)
 
if ( !is.null(args$help) || is.null(args$pdfinput) || is.null(args$noloutput) ) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}

library(DESeq2)
library(ggplot2)

rm(list=ls())

load("../../Figure2/Analysis.expression_normalize/deseq2.RData")

gene2type <- read.table('../../Figure2/Sup.Fig2.expression_boxplot/gene2type.xls',header=F)
mRNA <- gene2type[which(gene2type$V2=='Protein' & gene2type$V1 %in% rownames(assay(vsd))),]
lncRNA <- gene2type[which(gene2type$V2=='lncRNA' & gene2type$V1 %in% rownames(assay(vsd))),]
 
# PCA分析------------------------------------------------------mRNA

pcaData <- plotPCA(vsd, intgroup=c("Sample", "Tissue"), returnData=TRUE,plot3d = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Sample, shape=Tissue)) +
  geom_point(size=3) +
  scale_shape_manual(values = c(1:12)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

ggsave("mRNA.sample.pca.pdf")


# PCA分析------------------------------------------------------mRNA

pcaData <- plotPCA(assay(vsd)[lncRNA$V1,], intgroup=c("Sample", "Tissue"), returnData=TRUE, plot3d = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Sample, shape=Tissue)) +
  geom_point(size=3) +
  scale_shape_manual(values = c(1:12)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
pc <- prcomp(assay(vsd)[lncRNA$V1,])
print(pc)
ggsave("lncRNA.sample.pca.pdf")


# ------
install.packages('ggfortify')
library(ggfortify)
lncRNA_clean_data <- assay(vsd)[lncRNA$V1,]
samples <- sapply( strsplit(as.character(row.names(t(lncRNA_clean_data))), "-"), "[[", 1 )
metadata <- data.frame(samples,row.names = colnames(lncRNA_clean_data))
metadata$tissue <- vsd$Tissue

lncRNA_clean_data$group <- sapply( strsplit(as.character(row.names(lncRNA_clean_data)), "-"), "[[", 1 )
autoplot(object = prcomp(lncRNA_clean_data), data=lncRNA_clean_data, label=T, color=group)

mRNA_clean_data <- assay(vsd)[mRNA$V1,]
#lncRNA_clean_data$group <- sapply( strsplit(as.character(row.names(lncRNA_clean_data)), "-"), "[[", 1 )
autoplot(object = prcomp(t(assay(vsd)[lncRNA$V1,])), label=T, color=group)

# 还有一种方法 ---------------------------------
BiocManager::install('PCAtools')
library(PCAtools)
p <- pca(lncRNA_clean_data, metadata = metadata, removeVar = 0.1)

biplot(p,
       colby = 'samples', 
       lab=p$yvars,
       hline = 0, vline = c(-25, 0, 25),
       legendPosition = 'top', legendLabSize = 10, legendIconSize = 8.0)

ggsave('lncRNA.pca.pdf')

p <- pca(mRNA_clean_data, metadata = metadata, removeVar = 0.1)

biplot(p,
       colby = 'samples', 
       shape = 'tissue',
       lab=p$yvars,
       hline = 0, vline = c(-25, 0, 25),
       legendPosition = 'top', legendLabSize = 10, legendIconSize = 8.0)

ggsave('mRNA.pca.pdf')
