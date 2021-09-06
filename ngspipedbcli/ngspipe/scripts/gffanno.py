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

with open(sys.argv[3], 'r') as f:
    for line in f:
        line = line.rstrip()
        if line.startswith('#'):
            print(line)
            continue
        items = line.split('\t')
        feature_type, feature_attr = [items[2], items[8]]
        if feature_type == 'mRNA':
            matObj = re.match('ID=([^;]+)', feature_attr)
            if matObj:
                mRNA_new_id = matObj.group(1)
                anno = []
                for database in anno_df.columns:
                    db_anno = database + ':' + '{' + anno_df.loc[mRNA_new_id, database].replace(';',',') + '}'
                    anno.append(db_anno)
                line = line + ';' + 'Note=' + '+'.join(anno)
                print(line)
            else:
                sys.exit()
        else:
            print(line)


