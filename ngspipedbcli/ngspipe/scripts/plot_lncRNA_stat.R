# plot lncRNA

library(ggplot2)

len <- read.csv("/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/plot/transcript_length.csv",sep=",",header = FALSE)
exp <- read.csv("/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/assembly_final/gene.norm.tsv",sep="\t")
num <- read.csv("/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/plot/exon_num.csv",sep=",",header = FALSE)
classify <- read.csv("/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/cuffcompare_feature_info/cuffcompare_feature_info.tsv",sep="\t")

protein_transcript <- classify[which(classify$transcript_biotype == 'mRNA'), 'X']
lncRNA_transcript <- classify[which(classify$transcript_biotype %in% c('LincRNA','Processed transcript', 'Sense intronic','Antisense RNAs')), 'X']
TUCP <- classify[which(classify$transcript_biotype == 'TUCP'), 'X']

ggplot(len[protein_transcript,], aes(y=V2,x=1))+geom_boxplot()
ggplot(len[lncRNA_transcript,], aes(y=V2,x=2))+geom_boxplot()