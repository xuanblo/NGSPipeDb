# SimpleSQLite
# pip install simplesqlite
from simplesqlite import SimpleSQLite
import pandas as pd
import sys
import sqlite3

# read annotation file
print('reading annotation file...')
anno_df = pd.read_csv(sys.argv[1], sep='\t', index_col = 0, header = 0)

new_list = ['GeneID'] + list(anno_df.columns)

anno_df = anno_df.reindex(columns=new_list)

anno_df.GeneID = anno_df.index

print('connect gff sqlite3 database...')
con = SimpleSQLite(sys.argv[2])

print('create new table to gff sqlite3 database...')
con.create_table_from_dataframe(anno_df, table_name='anno')

print('create index to sqlite3 database...')
con.create_index_list('anno', ['GeneID', 'Nr', 'Swissprot', 'KEGG', 'TrEMBL', 'Interpro'])

con.commit()
con.close()