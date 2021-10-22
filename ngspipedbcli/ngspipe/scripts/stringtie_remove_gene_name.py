import sys
import pandas as pd

exp_matrix_with_gene_name_df = pd.read_csv(sys.argv[1], header=0, index_col=None)

exp_matrix_with_gene_name_df.iloc[:,0] = [i.split('|')[0] for i in exp_matrix_with_gene_name_df.iloc[:,0]]

exp_matrix_with_gene_name_df.to_csv(sys.argv[2], index=None)