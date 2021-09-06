import sys
import os
import pandas as pd

htseq = sys.argv[1:-1]

resultfile = sys.argv[-1]

df_arr = []
i = 0

for htseq_count in htseq:
    print(i, htseq_count)
    i = i + 1
    name = os.path.basename(htseq_count).split('.')[0]
    df_sample = pd.read_csv(htseq_count, index_col=0, header=None, sep = '\t')
    df_sample.columns = [name]
    bool_remove = df_sample.index.str.contains('__')
    df_sample = df_sample[~bool_remove]
    df_arr.append(df_sample)

df_all = pd.concat(df_arr, axis=1, sort=True)


df_all.to_csv(resultfile, sep='\t')