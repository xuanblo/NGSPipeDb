# 全基因功能分析
# 对所有unigenes进行功能注释和富集分析
# GO富集基因列表
# 柱状图
# KEGG富集基因列表
# 散点图
# 通路图
# 富集分析

# 使用clusterProfiler包对非模式生物进行GO和KEGG的富集分析
library(clusterProfiler)
# topGO

# download go term's name
# wget http://purl.obolibrary.org/obo/go/go-basic.obo

# download k's name
# https://www.kegg.jp/kegg/rest/keggapi.html
# http://rest.kegg.jp/find/ko/K15379
# http://rest.kegg.jp/get/ko:K00500

# 读取GO注释文件
go <- read.table("../result/eggnog/go_addname.txt",sep="\t",quote="")
bp = subset(go, go$V4=="biological_process")
mf = subset(go, go$V4=="molecular_function")

# 读取KEGG注释文件
kegg <- read.table("../result/eggnog/kegg_pathway_addname.txt",sep="\t",quote="")

# 读取差异基因列表
gene <- read.table("../result/diff/gene.tsv.0h_vs_11h.DESeq2.DE_results")
diff_gene = subset(gene,abs(log2FoldChange)>2&padj<0.05)
diffGeneList = row.names(diff_gene)

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

# 设置KEGG的TERM2GENE和设置TERM2NAME
kegg_term2gene <- data.frame(kegg$V4,kegg$V1)
kegg_term2name <- data.frame(kegg$V4,kegg$V5)
names(kegg_term2gene) <- c("ko_term","gene")
names(kegg_term2name) <- c("ko_term","name")

# 查看数据
head(bp_term2gene)
head(bp_term2name)
head(kegg_term2gene)
head(kegg_term2name)
gene <- as.vector(rownames(gene))
head(gene)

# 使用enrichr函数进行GO BP富集分析 (超几何分布)
bp_enrich <- enricher(gene=diffGeneList, pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                      TERM2GENE = bp_term2gene, TERM2NAME = bp_term2name, maxGSSize=NA)
write.csv(as.data.frame(bp_enrich),"../result/eggnog/BP_enrichment.csv",row.names = F)
# 使用enrichr函数进行GO MF富集分析 (超几何分布)
mf_enrich <- enricher(gene=diffGeneList, pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                      TERM2GENE = mf_term2gene, TERM2NAME = mf_term2name, maxGSSize=NA)
write.csv(as.data.frame(mf_enrich),"../result/eggnog/MF_enrichment.csv",row.names = F)

head(as.data.frame(bp_enrich))

# 合并mf和bp
df_bp = as.data.frame(bp_enrich)
df_mf = as.data.frame(mf_enrich)
df_bp$type="biological_process"
df_mf$type="molecular_function"
df_bp_mf = merge(df_bp, df_mf, all=TRUE)
df_bp_mf_sorted = df_bp_mf[sort(df_bp_mf$p.adjust,index.return=TRUE, decreasing = FALSE)$ix,]
df_bp_mf_finnal = df_bp_mf_sorted[0:20,]
write.csv(df_bp_mf_sorted,"../result/eggnog/BP_MF_enrichment.csv",row.names = F)

# 对富集结果进行可视化
library(ggplot2)
p = ggplot(df_bp_mf_finnal, aes(x=Description, y=Count, fill=type)) + geom_bar(stat='identity') + coord_flip()
p + labs(title = "The Most Enriched GO Terms", x = "", y = "Number of genes") + theme(legend.position="top")
ggsave("../result/eggnog/go_enrich.pdf")

# 对富集结果进行可视化
# 条形图
barplot(bp_enrich, showCategory=8)
barplot(mf_enrich, showCategory=8)

# 点图
dotplot(bp_enrich, x='Count')

## categorySize can be scaled by 'pvalue' or 'geneNum'

# 网络图
cnetplot(bp_enrich, categorySize="pvalue", foldChange=diffGeneList)
#emapplot(go_enrich)

# 使用enrichr函数进行KEGG的富集分析
# gene 待富集的基因list
kegg_enrich <- enricher(gene=diffGeneList, pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                        TERM2GENE = kegg_term2gene, TERM2NAME = kegg_term2name)
write.csv(as.data.frame(kegg_enrich),"../result/eggnog/KEGG_enrichment.csv",row.names = F)

head(as.data.frame(kegg_enrich))

# 对富集结果进行可视化
barplot(kegg_enrich, showCategory=8)
dotplot(kegg_enrich)
## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(kegg_enrich, categorySize="pvalue", foldChange=diffGeneList)
#emapplot(kegg_enrich)

# 使用pathview对感兴趣的基因画图
library(pathview)
pathview(gene.data = c("K00128", "K13953"), pathway.id = '00071', species = "ko", kegg.native = T)