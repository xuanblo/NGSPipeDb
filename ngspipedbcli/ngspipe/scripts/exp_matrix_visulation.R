# DESeq2

library(DESeq2)
library(ggplot2)
library(edgeR)

# https://www.jianshu.com/p/f4b618354dc2
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

# 常规表达矩阵，log2转换后或
# Deseq2的varianceStabilizingTransformation转换的数据
# 如果有批次效应，需要事先移除，可使用removeBatchEffect
# 如果有系统偏移(可用boxplot查看基因表达分布是否一致)，
# 需要quantile normalization
##导入数据##

# counts
dataExpr_norm <- read.table("/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/assembly_final/gene.norm.tsv", sep='\t', row.names=1, header=T, 
                       quote="", comment="", check.names=F)

dataExpr_orgin <- read.table("/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/assembly_final/gene.tsv", sep='\t', row.names=1, header=T, 
                             quote="", comment="", check.names=F)

# 读样本信息文件
coldata <- read.csv("/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/assembly_final/condition_keep.tsv", row.names=1,sep="\t")
coldata <- coldata[,c("condition","type")]

all(rownames(coldata) %in% colnames(countData))
countData <- countData[,rownames(coldata)]
all(rownames(coldata) == colnames(countData))

boxplot(log2(as.matrix(dataExpr_norm)+1), names=colnames(dataExpr), horizontal=TRUE)
boxplot(log2(as.matrix(dataExpr_orgin)+1), names=colnames(dataExpr_orgin), horizontal=TRUE)

# convert wide to long
library(tidyr)
dataExpr_orgin_df <- as.data.frame(log2(as.matrix(dataExpr_orgin)+1))
plotDat_orgin <- gather(dataExpr_orgin_df, "x", "y")
ggplot(data=plotDat_orgin, aes(x,y,color='red')) + geom_violin() + coord_flip()

dataExpr_norm_df <- as.data.frame(log2(as.matrix(dataExpr_norm)+1))
plotDat_norm <- gather(dataExpr_norm_df, "x", "y")
ggplot(data=plotDat_norm, aes(x,y,color='bule')) + geom_violin() + coord_flip()

# 过滤低表达量的基因
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 用rlog对count进行log转化（后续可视化分析）
# 方差稳定转换
rld <- rlog(dds, blind = FALSE) # 大于30样本比较慢
head(assay(rld), 3)
#vsd <- vst(dds, blind = FALSE) # 比较快
#head(assay(vsd), 3)

# 基因表达量方差过滤
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), dim(rld)[1]*0.25)

# 评估样品之间的总体相似性(距离)
#sampleDists <- dist(t(assay(rld)))
#sampleDists <- dist(t(norcounts))
#sampleDists <- dist(t(assay(vsd)))
sampleDists <- dist(t(assay(rld[topVarGenes,])))
#library("pheatmap")
#library("RColorBrewer")
#install.packages("hexbin")
sampleDistMatrix <- as.matrix( sampleDists )
# Blues
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
# filename = "/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/sample.heatmap.pdf"

# 评估样本之间相关性
cormat <- round(cor(assay(rld[topVarGenes,])),2)
#cormat <- round(cor(norcounts),2)
library(reshape2)

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
upper_tri
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# PCA分析

dev.new()
pcaData <- plotPCA(rld, intgroup=c("condition", "tissue"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=tissue)) +
  geom_point(size=3) +
  scale_shape_manual(values = c(1:12)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

ggsave("/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/sample.pca.pdf")
dev.off()

# ------------------

topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("condition"))

library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("condition"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = condition, y = count, color = condition)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
ggplot(geneCounts, aes(x = condition, y = count, color = condition, group = condition)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()

# 差异表达基因聚类，或者方差最大的基因聚类
library("genefilter")
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 20)
mat  <- assay(rld)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
pheatmap(mat, filename = "/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/varianceTop20.heatmap.pdf")


# tpm

tpm_orgin = read.table("/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/assembly_final/gene_tpm.txt", sep='\t', row.names=1, header=T, 
                       quote="", comment="", check.names=F)
boxplot(log2(as.matrix(tpm_orgin)+1), horizontal=TRUE)
# upper-quartile normalization
tpm_norm <- normalize.quantiles(tpm_orgin, copy = TRUE)
uq <- calcNormFactors(tpm_orgin, method="upperquartile")

# edge
lcpm <- cpm(dataExpr_orgin, log=TRUE)
keep.exprs <- rowSums(lcpm>1)>=3
lfcpm <- lcpm[keep.exprs,]
dim(lfcpm)

library(RColorBrewer)
nsamples <- ncol(lcpm)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

plot(density(lfcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lfcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")