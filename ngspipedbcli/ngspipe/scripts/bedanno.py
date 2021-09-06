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

#anno_df.to_csv(sys.argv[3], sep='\t')

with open(sys.argv[3], 'r') as f:
    for line in f:
        line = line.rstrip()
        items = line.split('\t')
        chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts = items
        score = "1000" # higher numbers = darker gray
        items[4] = score
        itemRgb = "128,0,128" # color
        items[8] = itemRgb
        blockSizes = blockSizes.rstrip(',')
        items[10] = blockSizes
        blockStarts = blockStarts.rstrip(',')
        items[11] = blockStarts
        geneName = name.split('.')[0]
        transcriptType = 'mRNA'
        anno = []
        for database in anno_df.columns:
            db_anno = database + ':' + '{' + anno_df.loc[name, database].replace(';',',').replace(' ', '_') + '}'
            #db_anno = database + ':' + '{' + anno_df.loc[name, database].replace(';',',') + '}'
            anno.append(db_anno)
        Annotation = 'Note=' + '+'.join(anno)
        items = items + [geneName, transcriptType, Annotation]
        line = '\t'.join(items)
        print(line)


