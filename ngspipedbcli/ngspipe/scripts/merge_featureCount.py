import sys
import os
import pandas as pd

featurecounts = sys.argv[1:-1]

resultfile = sys.argv[-1]


df_arr = []

j = 0

for featurecount in featurecounts:
    print(j, featurecount)
    j = j + 1
    name = os.path.basename(featurecount).split('.')[0]
    df_sample = pd.read_csv(featurecount, index_col=0, header=None, skiprows = 2, sep = '\t')
    df_sample.columns = ['a','b', 'c', 'd', 'e', name]
    df_arr.append(df_sample.loc[:,[name]])

df_all = pd.concat(df_arr, axis=1, sort=True)

df_all.to_csv(resultfile, sep='\t')