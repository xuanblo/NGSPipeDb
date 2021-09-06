#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 18:10:49 2021

@author: zhangxuan
"""

# test bioinfokit

from bioinfokit.analys import norm
from bioinfokit import analys, visuz

import pandas as pd
import re
import sys
import matplotlib.pyplot as plt

def count_bp(df):
    """Given a DataFrame with the exon coordinates from Gencode for a single
    gene, return the total number of coding bases in that gene.
    Example:
        >>> import numpy as np
        >>> n = 3
        >>> r = lambda x: np.random.sample(x) * 10
        >>> d = pd.DataFrame([np.sort([a,b]) for a,b in zip(r(n), r(n))], columns=['start','end']).astype(int)
        >>> d
           start  end
        0      6    9
        1      3    4
        2      4    9
        >>> count_bp(d)
        7
    Here is a visual representation of the 3 exons and the way they are added:
          123456789  Length
        0      ----       4
        1   --            2
        2    ------       6
            =======       7
    """
    start = df.start.min()
    end = df.end.max()
    bp = [False] * (end - start + 1)
    for i in range(df.shape[0]):
        s = df.iloc[i]['start'] - start
        e = df.iloc[i]['end'] - start + 1
        bp[s:e] = [True] * (e - s)
    return sum(bp)

anno = "../../../../testdata/GRCm38.83.chr19.gtf"
if anno.endswith('gtf'):
    GENE_ID = re.compile('gene_id "([^"]+)')
else:
    GENE_ID = re.compile('ID=([^;]+)')

exon_pos_dict = dict()

with open(anno, 'r') as f:
    for line in f:
        items = line.split('\t')
        if items[2] == 'exon':
            try:
                r = GENE_ID.search(items[8])
                gene_id = r.group(1)
                exon_pos_dict.setdefault(gene_id, {})
                exon_pos_dict[gene_id].setdefault('start', [])
                exon_pos_dict[gene_id].setdefault('end', [])
                exon_pos_dict[gene_id]['start'].append(int(items[3]))
                exon_pos_dict[gene_id]['end'].append(int(items[4]))
            except:
                sys.stderr.write("read gtf error\n")
                sys.exit(-1)

gene_len_dict = dict()

for g in exon_pos_dict.keys():
    start = exon_pos_dict[g]['start']
    end = exon_pos_dict[g]['end']
    gene_pos_df = pd.DataFrame({'start':start, 'end':end})
    length = count_bp(gene_pos_df)
    gene_len_dict[g] = length

gene_len_df = pd.DataFrame.from_dict(gene_len_dict,orient="index",columns=['Length'])

df = pd.read_csv('gene.csv', header=0, index_col=0, sep=',')
df.index = [i.split('|')[0] for i in df.index]

data = pd.concat([df, gene_len_df], axis=1)

# now, normalize raw counts using TPM method
# gene length must be in bp
nm = norm()
nm.tpm(df=data, gl='Length')
# get TPM normalized dataframe
tpm_df = nm.tpm_norm
tpm_df.head(2)
# output

# heatmap with hierarchical clustering 
visuz.gene_exp.hmap(df=tpm_df.iloc[0:100,:], cmap='RdYlGn', dim=(6, 6), tickfont=(6, 4), show=True)
