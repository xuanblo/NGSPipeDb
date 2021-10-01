import pandas as pd
import os
import glob
import sys

salmon_dir = sys.argv[1]

parent, samples, filenames = next(os.walk(salmon_dir))
matrix = pd.DataFrame()

for sample in sorted(samples):
    file_path = os.path.join(salmon_dir, sample, "quant.sf")
    df = pd.read_csv(file_path,index_col=0,header=0,sep='\t')
    matrix.loc[:,sample] = df.TPM

write2path = sys.argv[2]
matrix.to_csv(write2path, sep=',')
