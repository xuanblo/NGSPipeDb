import os
import sys
import pandas as pd

conditionfile = sys.argv[1]
qc_outdir = sys.argv[2]

df = pd.read_csv(conditionfile, header=0, index_col=None)

for col in df.index:
    print(df.iloc[col, 1], df.iloc[col, 0], os.path.join(qc_outdir, df.iloc[col, 0], df.iloc[col, 0]+'.cleanR1.fq.gz'), os.path.join(qc_outdir, df.iloc[col, 0], df.iloc[col, 0]+'.cleanR2.fq.gz'), sep='\t')
