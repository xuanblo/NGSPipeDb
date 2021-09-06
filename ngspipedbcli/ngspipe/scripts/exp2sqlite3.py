# SimpleSQLite
# pip install simplesqlite

from simplesqlite import SimpleSQLite
import pandas as pd
import sys
import sqlite3
import os

# read annotation file
print('reading expression file...')
anno_df = pd.read_csv(sys.argv[1], sep='\t', index_col = 0, header = 0)

new_list = ['Gene_id'] + list(anno_df.columns)

anno_df = anno_df.reindex(columns=new_list)

anno_df.Gene_id = anno_df.index

print('connect sqlite3 database...')
if os.path.exists(sys.argv[2]):
    os.remove(sys.argv[2])
    print("file exists! we will detele it.")
con = SimpleSQLite(sys.argv[2])

print('create new table to sqlite3 database...')
#con.create_table_from_dataframe(anno_df, table_name='exp',add_primary_key_column=1)

con.create_table_from_dataframe(anno_df, table_name='exp',primary_key='Gene_id')


# Migration
print('create index to sqlite3 database...')
con.create_index_list('exp', ['Gene_id',])

con.commit()
con.close()
