# DESeq2
# https://www.cnblogs.com/chenpeng1024/p/9260803.html

rm(list=ls())

load("../../../Figures/Figure2/expression_normalize/deseq2.RData")

# Differential expression analysis，差异表达分析的输入是原始reads文件，不能处理
## 运行deseq2
# 对原始dds进行标准化
dds <- DESeq(dds)
## 提取结果

diff_result_dir = "diff123/"

sample_comb <- combn(unique(samps$condition),2)

for (i in 1:dim(sample_comb)[2]){
  sample1 <- as.character(sample_comb[, i][1])
  sample2 <- as.character(sample_comb[, i][2])
  
  cat(sample1, "vs", sample2, ";")
  
  res <- results(dds, contrast=c("condition",sample2,sample1))
  
  # foldchange = log2(sample1 / sample2)
  # 要把表达量添加进去
  # https://blog.csdn.net/qazplm12_3/article/details/81221329?utm_source=blogxgwz1
  
  # 把表达量添加到结果中
  # 获得第一组数据均值
  base1 <- counts(dds, normalized=TRUE)[, colData(dds)$condition == sample1]
  if (is.vector(base1)){
    baseMean1 <- as.data.frame(base1)
  } else {
    baseMean1 <- as.data.frame(rowMeans(base1)) # rowMedians
  }
  colnames(baseMean1) <- sample1
  #head(baseMean1)
  
  # 获得第二组数据均值
  base2 <- counts(dds, normalized=TRUE)[, colData(dds)$condition == sample2]
  if (is.vector(base2)){
    baseMean2 <- as.data.frame(base2)
  } else {
    baseMean2 <- as.data.frame(rowMeans(base2)) # rowMedians
  }
  colnames(baseMean2) <- sample2
  #head(baseMean2)

  # 结果组合
  res <- cbind(baseMean1, baseMean2, as.data.frame(res))
  head(res)
  
  # 增加ID信息
  res <- cbind(ID=rownames(res), as.data.frame(res))
  res$baseMean <- rowMeans(cbind(base1, base2))
  
  # 校正后p-value为NA的复制为1
  res$padj[is.na(res$padj)] <- 1

  # 按pvalue排序, 把差异大的基因放前面
  res <- res[order(res$pvalue),]
  #head(res)
  
  # the base mean is the mean of normalized counts of all samples, normalizing for sequencing depth.
  resSig_all <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1, select=c('ID', sample1, sample2, 'log2FoldChange', 'padj'))
  #resSig_all_sorted <- resSig_all[ order(abs(resSig_all$log2FoldChange), decreasing = TRUE), ]
  write.csv(resSig_all, file= paste(diff_result_dir,sample1,"_vs_",sample2,".all.csv", sep=''), quote=F)
  
  cat("all:",dim(resSig_all)[1], ";")
  
  resSig_up <- subset(res, padj < 0.05 & log2FoldChange > 1, select=c('ID', sample1, sample2, 'log2FoldChange', 'padj'))
  #resSig_up_sorted <- resSig_up[ order(resSig_up$log2FoldChange, decreasing = TRUE), ]
  write.csv(resSig_up, file= paste(diff_result_dir,sample1,"_vs_",sample2,".up.csv", sep=''), quote=F)
  
  cat("up:", dim(resSig_up)[1], ";")
  
  resSig_down <- subset(res, padj < 0.05 & log2FoldChange < -1, select=c('ID', sample1, sample2, 'log2FoldChange', 'padj'))
  #resSig_down_sorted <- resSig_down[ order(resSig_down$log2FoldChange), ]
  write.csv(resSig_down, file= paste(diff_result_dir,sample1,"_vs_",sample2,".down.csv", sep=''), quote=F)
  
  cat("down:", dim(resSig_down)[1], '\n')
  
}

