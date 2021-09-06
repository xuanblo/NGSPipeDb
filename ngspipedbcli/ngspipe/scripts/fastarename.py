import sys
import pandas as pd
import re

anno_df = pd.read_excel(sys.argv[1], index_col=0, header=0)
anno_df = anno_df.fillna('No result')

mRNA_old2new_df = pd.read_csv(sys.argv[2], header=0, index_col=0, sep=',')

anno_new_index = []

for mRNA_old_id in anno_df.index:
    anno_new_index.append(mRNA_old2new_df.loc[mRNA_old_id, 'newid'])

anno_df.index = anno_new_index

anno_df.to_csv(sys.argv[3], sep='\t')

with open(sys.argv[4], 'r') as f:
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            matObj = re.match('>(xingyouteng_\d+).*', line)
            if matObj:
                mRNA_old_id = matObj.group(1)
                mRNA_new_id = mRNA_old2new_df.loc[mRNA_old_id, 'newid']
                line = line.replace(mRNA_old_id, mRNA_new_id)
                anno = []
                for database in anno_df.columns:
                    db_anno = database + ':' + '{' + anno_df.loc[mRNA_new_id, database] + '}'
                    anno.append(db_anno)
                line = line + '  ' + 'Note=' + ';'.join(anno)
                print(line)
        else:
            print(line)


