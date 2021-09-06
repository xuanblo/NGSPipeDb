#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 16:20:44 2019

@author: zhangxuan
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

length = pd.read_csv("/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/plot/transcript_length.csv",header = None, index_col=0)
exp = pd.read_csv("/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/assembly_final/gene.norm.tsv",sep="\t",header = 0, index_col=0)
num = pd.read_csv("/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/plot/exon_num.csv",sep=",",header = None, index_col=0)
classify = pd.read_csv("/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/cuffcompare_feature_info/cuffcompare_feature_info.tsv",sep="\t",header=0, index_col=0)

protein_transcript = classify[classify.transcript_biotype == 'mRNA'].index
lncRNA_transcript = classify[classify.transcript_biotype.isin(['LincRNA','Processed transcript', 'Sense intronic','Antisense RNAs'])].index
TUCP = classify[classify.transcript_biotype == 'TUCP'].index

# length

protein_transcript_length = list(np.log10(length.loc[protein_transcript,1] + 1))
lncRNA_transcript_length = list(np.log10(length.loc[lncRNA_transcript,1] + 1))
TUCP_transcript_length = list(np.log10(length.loc[TUCP,1] + 1))
sns.boxplot([protein_transcript_length,lncRNA_transcript_length,TUCP_transcript_length])

plt.figure()
plt.boxplot([protein_transcript_length,lncRNA_transcript_length], labels=['protein','lncRNA'])
plt.ylabel('transcript length log10+1')
plt.show()


# exon number

protein_transcript_exon_num = list(np.log2(num.loc[protein_transcript,1]))
lncRNA_transcript_exon_num = list(np.log2(num.loc[lncRNA_transcript,1]))
TUCP_transcript_exon_num = list(np.log2(num.loc[TUCP,1]))

plt.figure()
plt.boxplot([protein_transcript_exon_num,lncRNA_transcript_exon_num], labels=['protein','lncRNA'])
plt.ylabel('transcript exon number log2')
plt.show()


# expression

protein_transcript_exp_num = (exp>0.1).sum()
plt.bar(protein_transcript_exp_num.index, protein_transcript_exp_num.values)