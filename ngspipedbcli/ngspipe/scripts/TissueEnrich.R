args <- commandArgs( trailingOnly = TRUE )

input_exp_matrix = args[1]
output_ts_gene = args[2]

library(TissueEnrich)
# 默认FPKM过滤是多少

# Tissue-specific gene retrieval
expressionData<-read.table(input_exp_matrix,header=TRUE,row.names=1,sep='\t')

se<-SummarizedExperiment(assays = SimpleList(as.matrix(expressionData)),
                         rowData = row.names(expressionData),colData = colnames(expressionData))
output<-teGeneRetrieval(se)
head(assay(output))
df = data.frame(assay(output))
df_ts = subset(df,df$Group=="Tissue-Enriched"|df$Group=="Group-Enriched"|df$Group=="Tissue-Enhanced")

#The genes are divided into six groups based on their gene expression across the tissues. These groups are:
  
# - Not Expressed: Genes with an expression level less than 1 (TPM or FPKM) across all the tissues.
# - Tissue Enriched: Genes with an expression level greater than or equal to 1 (TPM or FPKM) 
    # that also have at least five-fold higher expression levels in a particular tissue compared to all other tissues.
# - Group Enriched: Genes with an expression level greater than or equal to 1 (TPM or FPKM) 
    # that also have at least five-fold higher expression levels in a group of 2-7 tissues compared to all other tissues, 
    # and that are not considered Tissue Enriched.
# - Tissue Enhanced: Genes with an expression level greater than or equal to 1 (TPM or FPKM) 
    # that also have at least five-fold higher expression levels in a particular tissue compared to the average levels 
    # in all other tissues, and that are not considered Tissue Enriched or Group Enriched.
# - Expressed in all: Genes with an expression level greater than or equal to 1 (TPM or FPKM) across all of the tissues 
    # that are not in any of the above 4 groups.
# - Mixed: Genes that are not assigned to any of the above 5 groups.
# Genes from the Tissue Enriched, Group Enriched, and Tissue Enhanced groups are classified as tissue-specific genes.

write.table(df_ts, output_ts_gene, row.names = FALSE)