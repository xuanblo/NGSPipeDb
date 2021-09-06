import sys
import re
import pandas as pd

species_prefix = 'Pvo'
chr = 'chr1'
gene_num = 0

def new_gene_id(gene_old_id='nbis_noL1id-mrna-14664'):
    gene_old_id = gene_old_id[17:]

gene_old2new = dict()
mRNA_old2new = dict()
g2t_num = dict()
t2exon_num = dict()
t2cds_num = dict()

with open(sys.argv[1], 'r') as f:
    for line in f:
        line = line.rstrip()
        if line.startswith('#'):
            print(line)
        else:
            items = line.split('\t')
            features_type = items[2]
            features_attr = items[8]
            features_loc = items[0]
            if features_loc.startswith('chr'):
                chr = features_loc.lstrip('chr')
            else:
                chr = 'X'
            items[1] = 'PVDB'
            items[5] = '.'
            line = '\t'.join(items)
            if features_type == 'gene':
                matObj = re.match('ID=(.*)', features_attr)
                if matObj:
                    gene_old_id = matObj.group(1)
                    gene_num = gene_num + 1
                    gene_id_tmp = str(gene_num)
                    gene_id_tmp = gene_id_tmp.zfill(6)
                    gene_new_id = species_prefix + chr + 'g' + gene_id_tmp
                    gene_old2new[gene_old_id] = gene_new_id
                    line = line.replace(gene_old_id, gene_new_id)
                    print(line)
            elif features_type == 'mRNA':
                matObj = re.match('ID=([^;]+);Parent=(.*)', features_attr)
                if matObj:
                    mRNA_old_id, gene_old_id = matObj.groups()
                    gene_new_id = gene_old2new[gene_old_id]
                    g2t_num.setdefault(gene_new_id, 0)
                    g2t_num[gene_new_id] += 1
                    mRNA_new_id = gene_new_id + '.{}'.format(g2t_num[gene_new_id]) # mRNA
                    mRNA_old2new[mRNA_old_id] = mRNA_new_id
                    line = line.replace(gene_old_id, gene_new_id)
                    line = line.replace(mRNA_old_id, mRNA_new_id)
                    print(line)
            elif features_type == 'exon':
                matObj = re.match('ID=([^;]+);Parent=(.*)', features_attr)
                if matObj:
                    exon_old_id, mRNA_old_id = matObj.groups()
                    mRNA_new_id = mRNA_old2new[mRNA_old_id]
                    t2exon_num.setdefault(mRNA_new_id, 0)
                    t2exon_num[mRNA_new_id] += 1
                    exon_new_id = mRNA_new_id + ':exon:{}'.format(t2exon_num[mRNA_new_id])
                    mRNA_old2new[mRNA_old_id] = mRNA_new_id
                    line = line.replace(exon_old_id, exon_new_id)
                    line = line.replace(mRNA_old_id, mRNA_new_id)
                    print(line)
            elif features_type == 'CDS':
                matObj = re.match('ID=([^;]+);Parent=(.*)', features_attr)
                if matObj:
                    cds_old_id, mRNA_old_id = matObj.groups()
                    mRNA_new_id = mRNA_old2new[mRNA_old_id]
                    t2cds_num.setdefault(mRNA_new_id, 0)
                    t2cds_num[mRNA_new_id] += 1
                    cds_new_id = mRNA_new_id + ':cds:{}'.format(t2cds_num[mRNA_new_id])
                    mRNA_old2new[mRNA_old_id] = mRNA_new_id
                    line = line.replace(cds_old_id, cds_new_id)
                    line = line.replace(mRNA_old_id, mRNA_new_id)
                    print(line)

df = pd.DataFrame.from_dict(mRNA_old2new, orient='index',columns=['newid'])
df.to_csv(sys.argv[2])



