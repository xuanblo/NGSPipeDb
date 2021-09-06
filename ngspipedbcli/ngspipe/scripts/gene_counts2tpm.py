import sys
import os
import re
import pandas as pd
import numpy as np

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

def read_counts2tpm(df, sample_name):
    """
    convert read counts to TPM (transcripts per million)
    :param df: a dataFrame contains the result coming from featureCounts
    :param sample_name: a list, all sample names, same as the result of featureCounts
    :return: TPM
    """
    result = df
    sample_reads = result.loc[:, sample_name].copy()
    gene_len = result.loc[:, ['Length']]
    rate = sample_reads.values / gene_len.values
    tpm = rate / np.sum(rate, axis=0).reshape(1, -1) * 1e6
    return pd.DataFrame(data=tpm, index=df.index, columns=sample_name)

# read gtf

GENE_ID = re.compile('gene_id "([^"]+)"')

exon_pos_dict = dict()

with open(sys.argv[1], 'r') as f:
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

# read featurecounts reasult
featureCounts = pd.read_csv(sys.argv[2], header=1, index_col=0, sep='\t')
featureCounts = featureCounts.iloc[:,5:]
featureCounts.columns = [os.path.basename(i).split('.')[0] for i in featureCounts.columns]

data = pd.concat([featureCounts, gene_len_df], axis=1)

df = read_counts2tpm(data, featureCounts.columns)

df.to_csv(sys.argv[3], sep='\t')