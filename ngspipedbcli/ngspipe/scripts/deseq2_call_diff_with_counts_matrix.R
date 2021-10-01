#!/usr/env/bin/Rscript
# _*_ coding: utf-8 _*_

#options(getopt.quiet = TRUE)
if (!require('getopt', character.only = TRUE)) {
  install.packages('getopt', dependencies = TRUE, quiet = TRUE)
  #suppressMessages(library('getopt', character.only = TRUE))
  library('getopt', character.only = TRUE)
}
command = matrix(c(
  'help'           , 'h', 0, "logical", 'help message',
  'countsMatrix'   , 'm', 1, "character", 'expression must be counts',
  'conditionFile'  , 'c', 1, "character",'sample group information',
  'resultDir'      , 'r', 1, "character",'sample compare result directory'
), byrow=TRUE, ncol=5)
args = getopt(command)

if ( !is.null(args$help) || is.null(args$countsMatrix) || is.null(args$conditionFile) || is.null(args$resultDir) ) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}

## 设置默认值
if ( is.null(args$countsMatrix) ) {
  args$countsMatrix = ""
  }

chooseCRANmirror(ind=16)
# install normal packages
my_packages_normal <- c("getopt", "ggplot2", "BiocManager")                                        # Specify your packages
not_installed_normal <- my_packages_normal[!(my_packages_normal %in% installed.packages()[ , "Package"])]    # Extract not installed packages
if(length(not_installed_normal)) install.packages(not_installed_normal)                               # Install not installed packages

# install bioconductor packages
my_packages_bioconductor <- c("GenomeInfoDb", "DESeq2")                                        # Specify your packages
not_installed_bioconductor <- my_packages_bioconductor[!(my_packages_bioconductor %in% installed.packages()[ , "Package"])]    # Extract not installed packages
if(length(not_installed_bioconductor)) install.packages(not_installed_bioconductor)                               # Install not installed packages

#suppressMessages(library('DESeq2'))
#suppressMessages(library('ggplot2'))
library('DESeq2')
#deseq2需要GenomeInfoDb, GenomeInfoDb又需要GenomeInfoDbData，conda在macos上安装不了GenomeInfoDbData
# /Users/zhangxuan/opt/anaconda3/envs/exp_r_env/lib/R/library/


# DESeq2
# https://www.cnblogs.com/chenpeng1024/p/9260803.html


# 读表达量矩阵
# args$countsMatrix="/Users/zhangxuan/Work/Current_work2020-6-21/databasetool/mouse_transcriptome_analysis/results/result/quantify/quantify_by_stringtie/gene.csv"
cts <- read.csv(args$countsMatrix, sep=",", row.names=1, check.names=FALSE)
colnames(cts) <- sub("\\.", "-", colnames(cts))

countData <- as.matrix(cts)
rownames(countData) <- rownames(cts)
countData <- countData[,sort(colnames(countData),decreasing = FALSE)]

# 读样本信息文件
# args$conditionFile="/Users/zhangxuan/Work/Current_work2020-6-21/databasetool/mouse_transcriptome_analysis/testdata/condition.xls"
coldata <- read.csv(args$conditionFile, row.names=1, sep=",", check.names=FALSE)
coldata <- coldata[,c("Sample","Tissue")]
coldata <- coldata[sort(rownames(coldata),decreasing = FALSE),]

# sort rows and columns
if(!all(rownames(coldata) %in% colnames(countData)))
{
  print("error: not all samples in exp matrix and condition are the same")
  q("no", 1, FALSE)
} else {
  if(!all(rownames(coldata) == colnames(countData)))
  {
    countData <- countData[, rownames(coldata)]
  }
}

# 创建deseq2对象
dds <- DESeqDataSetFromMatrix(countData = round(countData),
                              colData = coldata,
                              design = ~ Sample)
dds

# 过滤低表达量的基因
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Differential expression analysis，差异表达分析的输入是原始reads文件，不能处理
## 运行deseq2
# 对原始dds进行标准化
dds <- DESeq(dds)
## 提取结果

#args$resultDir = "/Users/zhangxuan/Work/Current_work2020-6-21/databasetool/mouse_transcriptome_analysis/results/result/diff_samples"
args$resultDir = paste(args$resultDir,'/',sep="")

# 输出表达量均一化
#dds <- estimateSizeFactors(dds)
#sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

write.csv(normalized_counts, paste(args$resultDir, "/normalized.counts.csv", sep=""), quote=F)

save(dds, coldata, file = paste(args$resultDir, "/deseq2.RData" ,sep=""))
#save.image(paste(args$resultDir, "/deseq2.RData" ,sep=""))

sample_comb <- combn(unique(coldata$Sample),2)

detail_sample_name = function (sample){
  paste(sample, 'median', sep='_')
}

for (i in 1:dim(sample_comb)[2]){
  sample1 <- as.character(sample_comb[, i][1])
  sample2 <- as.character(sample_comb[, i][2])
  
  cat(sample1, "vs", sample2, ";")
  
  res <- results(dds, alpha=0.05, contrast=c("Sample",sample2,sample1))
  
  # foldchange = log2(sample1 / sample2)
  # 要把表达量添加进去
  # https://blog.csdn.net/qazplm12_3/article/details/81221329?utm_source=blogxgwz1
  
  # 把表达量添加到结果中
  # 获得第一组数据均值
  base1 <- counts(dds, normalized=TRUE)[, colData(dds)$Sample == sample1]
  if (is.vector(base1)){
    baseMean1 <- as.data.frame(base1)
  } else {
    baseMean1 <- as.data.frame(rowMedians(base1)) # rowMedians or rowMeans
  }
  colnames(baseMean1) <- detail_sample_name(sample1)
  #head(baseMean1)
  
  # 获得第二组数据均值
  base2 <- counts(dds, normalized=TRUE)[, colData(dds)$Sample == sample2]
  if (is.vector(base2)){
    baseMean2 <- as.data.frame(base2)
  } else {
    baseMean2 <- as.data.frame(rowMedians(base2)) # rowMedians or rowMeans
  }
  colnames(baseMean2) <- detail_sample_name(sample2)
  #head(baseMean2)
  
  # 结果组合
  ## add gene expression value of all samples 
  res <- cbind(baseMean1, as.data.frame(base1), baseMean2, as.data.frame(base2), as.data.frame(res))
  #head(res)
  
  # 增加ID信息
  res <- cbind(ID=rownames(res), as.data.frame(res))
  res$baseMean <- rowMeans(cbind(base1, base2))
  
  # 校正后p-value为NA的复制为1
  res$padj[is.na(res$padj)] <- 1
  
  # 按pvalue排序, 把差异大的基因放前面
  res <- res[order(res$pvalue),]
  #head(res)
  
  # the base mean is the mean of normalized counts of all samples, normalizing for sequencing depth.
  #resSig_all <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1, select=c('ID', sample1, sample2, 'log2FoldChange', 'padj'))
  resSig_all <- subset(res, select=c('ID', detail_sample_name(sample1), colnames(base1), detail_sample_name(sample2), colnames(base2), 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'))
  #colnames(resSig_all) = c('ID', detail_sample_name(sample1), detail_sample_name(sample2), 'log2FoldChange', 'padj')
  write.table(resSig_all, file= paste(args$resultDir,sample2,"_vs_",sample1,".all.csv", sep=''), quote=F, row.names = FALSE, sep=',')
  
  cat("all:",dim(resSig_all)[1], ";")
  
  resSig_up <- subset(res, padj < 0.05 & log2FoldChange > 1, select=c('ID', detail_sample_name(sample1), colnames(base1), detail_sample_name(sample2), colnames(base2), 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'))
  #colnames(resSig_up) = c('ID', detail_sample_name(sample1), detail_sample_name(sample2), 'log2FoldChange', 'padj')
  write.table(resSig_up, file= paste(args$resultDir,sample2,"_vs_",sample1,".up.padj.csv", sep=''), quote=F, row.names = FALSE, sep=',')
  
  cat("up_padj<0.05_fc>2:", dim(resSig_up)[1], ";")
  
  resSig_down <- subset(res, padj < 0.05 & log2FoldChange < -1, select=c('ID', detail_sample_name(sample1), colnames(base1), detail_sample_name(sample2), colnames(base2), 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'))
  #colnames(resSig_down) = c('ID', detail_sample_name(sample1), detail_sample_name(sample2), 'log2FoldChange', 'padj')
  write.table(resSig_down, file= paste(args$resultDir,sample2,"_vs_",sample1,".down.padj.csv", sep=''), quote=F, row.names = FALSE, sep=',')
  
  cat("down_padj<0.05_fc>2:", dim(resSig_down)[1], '\n')
}


# write.csv(resdata, "all_des_output.csv", row.names=FALSE)

