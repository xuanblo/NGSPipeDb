import sys
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np

import matplotlib
matplotlib.use('Agg')

gene_exp_matrix_file = sys.argv[1]

out_pic_file = sys.argv[2]


exp_df = pd.read_csv(gene_exp_matrix_file, header=0, index_col=0, sep='\t')

exp_df = exp_df.loc[~(exp_df==0).all(axis=1), :]  # 删了它


plt.figure()
plt.boxplot(np.log2(np.array(exp_df.values)+1))
plt.xticks(range(len(exp_df.columns)),exp_df.columns)
plt.ylabel("log2 FPKM+1")
plt.savefig(out_pic_file, format='pdf', bbox_inches='tight')
#plt.show()
plt.close()