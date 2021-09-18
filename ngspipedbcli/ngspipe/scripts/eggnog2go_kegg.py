import sys
import re

gtf_file = sys.argv[1]
eggnog = sys.argv[2]
gofile = sys.argv[3]
kofile = sys.argv[4]
Kfile = sys.argv[5]

rna2gene = dict()
gene2go = dict()
gene2ko = dict()
gene2K = dict()

if gtf_file.endswith('gtf'):
    gtf_gene_pat = re.compile('gene_id "([^"]+)')
    gtf_transcript_pat = re.compile('gene_id "([^"]+).*transcript_id "([^"]+)')
    gtf_transcript_pat_rev = re.compile('transcript_id "([^"]+).*gene_id "([^"]+)')
else:
    gtf_gene_pat = re.compile('ID=([^;]+)')
    gtf_transcript_pat = re.compile('ID=([^;]+).*Parent=([^;]+)')
    gtf_transcript_pat_rev = re.compile('Parent=([^;]+).*ID=([^;]+)')

with open(gtf_file, 'r') as f:
    for line in f:
        line = line.rstrip()
        items = line.split('\t')
        matObj = gtf_transcript_pat.match(items[8])
        matObj_rev = gtf_transcript_pat_rev.match(items[8])
        #print(items[8])
        tid = ''
        gid = ''
        if matObj:
            if gtf_file.endswith('gtf'):
                gid = matObj.group(1)
                tid = matObj.group(2)
            else:
                gid = matObj.group(2)
                tid = matObj.group(1)
        elif matObj_rev:
            if gtf_file.endswith('gtf'):
                gid = matObj_rev.group(2)
                tid = matObj_rev.group(1)
            else:
                gid = matObj_rev.group(1)
                tid = matObj_rev.group(2)
        if tid:
            rna2gene[tid] = gid

with open(eggnog, 'r') as f:
    for line in f:
        if line.startswith("#"):
            continue
        else:
            items = line.split('\t')
            rna = items[0]
            if rna in rna2gene.keys():
                gene = rna2gene[rna]
                GO_terms = items[9]
                KEGG_KOs = items[12]
                KEGG_Ks = items[11]

                if GO_terms:
                    for term in GO_terms.split(','):
                        gene2go.setdefault(gene, {})
                        gene2go[gene][term] = ''
                if KEGG_KOs:
                    for ko in KEGG_KOs.split(','):
                        if 'ko' in ko:
                            gene2ko.setdefault(gene, {})
                            gene2ko[gene][ko] = ''
                if KEGG_Ks:
                    for K in KEGG_Ks.split(','):
                        if 'ko' in K:
                            gene2K.setdefault(gene, {})
                            gene2K[gene][K[3:]] = ''

with open(gofile, 'w') as f:
    for gene in gene2go.keys():
        for term in gene2go[gene].keys():
            f.write(gene + '\t' + term + '\n')

with open(kofile, 'w') as f:
    for gene in gene2ko.keys():
        for ko in gene2ko[gene].keys():
            f.write(gene + '\t' + ko + '\n')

with open(Kfile, 'w') as f:
    for gene in gene2K.keys():
        for ko in gene2K[gene].keys():
            f.write(gene + '\t' + ko + '\n')

