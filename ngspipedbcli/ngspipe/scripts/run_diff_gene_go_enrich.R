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

# go enrichment analysis
#BiocManager::install('clusterProfiler')
library(clusterProfiler)

# 读取GO注释文件
#go <- read.table("interpro.goslim.go.addname.txt",sep="\t",quote="")
go <- read.table(args$term_anno, sep="\t", quote="")
bp = subset(go, go$V4=="biological_process")
mf = subset(go, go$V4=="molecular_function")
cc = subset(go, go$V4=="cellular_component")

# 读取差异基因列表
gene <- read.table(args$difflist, sep=',', header = T)

#rownames(gene) <- data.frame(matrix(unlist(strsplit(as.character(gene$ID),'\\|')),ncol = 2, byrow = TRUE))$X1
#diff_gene = subset(gene,abs(log2FoldChange)>2&padj<0.05)
#diffGeneList = row.names(gene)
diffGeneList <- intersect(row.names(gene), row.names(go))

# 设置GO BP的TERM2GENE和设置TERM2NAME
bp_term2gene <- data.frame(bp$V2,bp$V1)
bp_term2name <- data.frame(bp$V2,bp$V3)
names(bp_term2gene) <- c("bp_term","gene")
names(bp_term2name) <- c("bp_term","name")

# 设置GO MF的TERM2GENE和设置TERM2NAME
mf_term2gene <- data.frame(mf$V2,mf$V1)
mf_term2name <- data.frame(mf$V2,mf$V3)
names(mf_term2gene) <- c("mf_term","gene")
names(mf_term2name) <- c("mf_term","name")

# 设置GO CC的TERM2GENE和设置TERM2NAME
cc_term2gene <- data.frame(cc$V2,cc$V1)
cc_term2name <- data.frame(cc$V2,cc$V3)
names(cc_term2gene) <- c("cc_term","gene")
names(cc_term2name) <- c("cc_term","name")

# 查看数据
#head(bp_term2gene)
#head(bp_term2name)


# 使用enrichr函数进行GO BP富集分析 (超几何分布)
bp_enrich <- enricher(gene=diffGeneList, pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                      TERM2GENE = bp_term2gene, TERM2NAME = bp_term2name, maxGSSize=NA)
#write.csv(as.data.frame(bp_enrich),"BP_enrichment.csv",row.names = F)
# 使用enrichr函数进行GO MF富集分析 (超几何分布)
mf_enrich <- enricher(gene=diffGeneList, pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                      TERM2GENE = mf_term2gene, TERM2NAME = mf_term2name, maxGSSize=NA)
#write.csv(as.data.frame(mf_enrich),"MF_enrichment.csv",row.names = F)
cc_enrich <- enricher(gene=diffGeneList, pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                      TERM2GENE = cc_term2gene, TERM2NAME = cc_term2name, maxGSSize=NA)

#head(as.data.frame(bp_enrich))

# 合并mf和bp
df_bp = as.data.frame(bp_enrich)
df_mf = as.data.frame(mf_enrich)
df_cc = as.data.frame(cc_enrich)

if (nrow(df_bp)!=0){df_bp$type="biological_process"}
if (nrow(df_mf)!=0){df_mf$type="molecular_function"}
if (nrow(df_cc)!=0){df_cc$type="cellular_component"}

if (nrow(df_bp)!=0){
    df_bp_mf = df_bp
    if (nrow(df_mf)!=0){
        if (nrow(df_cc)!=0){
            df_bp_mf = merge(c(df_bp, df_mf, df_cc), all=TRUE)
        }
        else{
            df_bp_mf = merge(c(df_bp, df_mf), all=TRUE)
        }
    }
    else{
        if (nrow(df_cc)!=0){
            df_bp_mf = merge(c(df_bp, df_cc), all=TRUE)
        }
    }
}else{
    if (nrow(df_mf)!=0){
        df_bp_mf = df_mf
        if (nrow(df_cc)!=0){
            df_bp_mf = merge(c(df_mf, df_cc), all=TRUE)
        }
    }
    else{
        if (nrow(df_cc)!=0){
            df_bp_mf = df_cc
        }
        else{
            df_bp_mf = df_bp
        }
    }
}

#df_bp_mf = merge(c(df_bp, df_mf, df_cc), all=TRUE)
df_bp_mf_sorted = df_bp_mf[sort(df_bp_mf$p.adjust,index.return=TRUE, decreasing = FALSE)$ix,]
df_bp_mf_finnal = df_bp_mf_sorted[,]
write.csv(df_bp_mf_sorted,args$output,row.names = F)

# 对富集结果进行可视化
library(ggplot2)
#p = ggplot(df_bp_mf_finnal, aes(x=Description, y=Count, fill=type)) + geom_bar(stat='identity') + coord_flip()
p = ggplot(df_bp_mf_finnal, aes(x=Description, y=Count)) + geom_bar(stat='identity') + coord_flip()
p + labs(title = "The Most Enriched GO Terms", x = "", y = "Number of genes") + theme(legend.position="top")
ggsave(paste(args$output, '.pdf', sep=''))
