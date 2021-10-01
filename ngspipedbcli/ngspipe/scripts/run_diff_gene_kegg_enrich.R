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
  'term_anno'   , 't', 1, "character", 'go or kegg annotation with description inside',
  'output'  , 'o', 1, "character",'enriched table file',
  'difflist', 'd', 1, "character",'differential expression gene list'
), byrow=TRUE, ncol=5)
args = getopt(command)

if ( !is.null(args$help) || is.null(args$term_anno) || is.null(args$output) || is.null(args$difflist)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}

# kegg enrichment analysis
#BiocManager::install('clusterProfiler')
library(clusterProfiler)

# 读取KEGG注释文件
#go <- read.table("interpro.goslim.go.addname.txt",sep="\t",quote="")
kegg <- read.table(args$term_anno, sep="\t", quote="")

# 读取差异基因列表
gene <- read.table(args$difflist, sep=',', header = T)
#rownames(gene) <- data.frame(matrix(unlist(strsplit(as.character(gene$ID),'\\|')),ncol = 2, byrow = TRUE))$X1
#diff_gene = subset(gene,abs(log2FoldChange)>2&padj<0.05)
#diffGeneList <- row.names(gene)

# 设置GO BP的TERM2GENE和设置TERM2NAME
kegg_term2gene <- data.frame(kegg$V2,kegg$V1)
kegg_term2name <- data.frame(kegg$V2,kegg$V3)
names(kegg_term2gene) <- c("kegg_term","gene")
names(kegg_term2name) <- c("kegg_term","name")

# 查看数据
#head(kegg_term2gene)
#head(kegg_term2name)
diffGeneList <- intersect(row.names(gene), row.names(kegg))

# 使用enrichr函数进行kegg富集分析 (超几何分布)
kegg_enrich <- enricher(gene=diffGeneList, pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                      TERM2GENE = kegg_term2gene, TERM2NAME = kegg_term2name, maxGSSize=NA)
#write.csv(as.data.frame(bp_enrich),"BP_enrichment.csv",row.names = F)

#head(as.data.frame(kegg_enrich))


df_bp_mf = kegg_enrich
df_bp_mf_sorted = df_bp_mf[sort(df_bp_mf$p.adjust,index.return=TRUE, decreasing = FALSE)$ix,]
df_bp_mf_finnal = df_bp_mf_sorted[,]
write.csv(df_bp_mf_sorted,args$output,row.names = F)

# 对富集结果进行可视化
library(ggplot2)
p = ggplot(df_bp_mf_finnal, aes(x=Description, y=Count)) + geom_bar(stat='identity') + coord_flip()
p + labs(title = "The Most Enriched KO Terms", x = "", y = "Number of genes") + theme(legend.position="top")
ggsave(paste(args$output, '.pdf', sep=''))
