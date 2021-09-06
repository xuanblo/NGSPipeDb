# DESeq2

library(DESeq2)
library(ggplot2)

# 读表达量矩阵
cts <- read.csv("/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/assembly_final/gene.tsv",sep="\t",row.names="Geneid")
colnames(cts) <- sub("\\.", "-", colnames(cts))

countData <- as.matrix(cts)
rownames(countData) <- rownames(cts)

# 读样本信息文件
coldata <- read.csv("/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/assembly_final/condition.tsv", row.names=1,sep="\t")
coldata <- coldata[,c("condition","type")]

all(rownames(coldata) %in% colnames(countData))

all(rownames(coldata) == colnames(countData))

#countData <- countData[,rownames(coldata)]

# 创建deseq2对象
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = coldata,
                              design = ~ condition)
dds

## 获得normalize的count
dds <- estimateSizeFactors(dds)
norcounts <- counts(dds, normalized=T)
write.table(norcounts, file="/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/assembly_final/gene.norm.tsv",
            sep='\t', quote = FALSE)