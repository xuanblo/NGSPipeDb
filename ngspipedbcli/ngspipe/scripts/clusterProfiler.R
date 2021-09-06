# 安装
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

# 使用clusterProfiler包对非模式生物进行GO和KEGG的富集分析
library(clusterProfiler)
getwd()
# 读取GO注释文件
go <- read.table("test_go_annot.txt",sep="\t",quote="")
# 读取KEGG注释文件
kegg <- read.table("test_kegg_annot.txt",sep="\t",quote="")
# 读取差异基因列表
gene <- read.table("test_gene.txt")
head(go)
head(kegg)
head(gene)
# 设置TERM2GENE
go_term2gene <- data.frame(go$V1,go$V4)
# 设置TERM2NAME
go_term2name <- data.frame(go$V1,go$V2)
names(go_term2gene) <- c("go_term","gene")
names(go_term2name) <- c("go_term","name")

kegg_term2gene <- data.frame(kegg$V2,kegg$V1)
kegg_term2name <- data.frame(kegg$V2,kegg$V3)
names(kegg_term2gene) <- c("ko_term","gene")
names(kegg_term2name) <- c("ko_term","name")
# 查看数据
head(go_term2gene)
head(go_term2name)
head(kegg_term2gene)
head(kegg_term2name)
gene <- as.vector(gene$V1)
head(gene)

# 使用enrichr函数进行GO富集分析
go_enrich <- enricher(gene=gene,pvalueCutoff = 0.05,pAdjustMethod = "BH",TERM2GENE = go_term2gene,TERM2NAME = go_term2name)
head(as.data.frame(go_enrich))
write.csv(as.data.frame(go_enrich),"GO_enrichment.csv",row.names = F)
# 对富集结果进行可视化
barplot(go_enrich, showCategory=8)
dotplot(go_enrich)
## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(go_enrich, categorySize="pvalue")
#emapplot(go_enrich)

# 使用enrichr函数进行KEGG的富集分析
kegg_enrich <- enricher(gene=gene,pvalueCutoff = 0.05,pAdjustMethod = "BH",TERM2GENE = kegg_term2gene,TERM2NAME = kegg_term2name)
head(as.data.frame(kegg_enrich))
write.csv(as.data.frame(kegg_enrich),"KEGG_enrichment.csv",row.names = F)
barplot(kegg_enrich, showCategory=8)
dotplot(kegg_enrich)
## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(kegg_enrich, categorySize="pvalue")
#emapplot(kegg_enrich)