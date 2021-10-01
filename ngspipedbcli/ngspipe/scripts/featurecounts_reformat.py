import sys
import os
import pandas as pd

featurecounts_result = sys.argv[1]

df = pd.read_csv(featurecounts_result, header=1, index_col=0, sep='\t')

df = df.iloc[:, 5:]

df.columns = [os.path.basename(os.path.dirname(i)) for i in df.columns]

df.to_csv(sys.argv[2])