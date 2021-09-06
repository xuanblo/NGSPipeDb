# make network object and load DataFrame, df
import sys
import pandas as pd
from clustergrammer import Network
df = pd.read_csv(sys.argv[1], header=True, index_col=0, sep='\t')
net = Network()
net.load_df(df)

# Z-score normalize the rows
net.normalize(axis='row', norm_type='zscore', keep_orig=True)

# filter for the top 100 columns based on their absolute value sum
net.filter_N_top('col', 100, 'sum')

# cluster using default parameters
net.cluster()

# save visualization JSON to file for use by front end
net.write_json_to_file('viz', sys.argv[2])