#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 22:13:16 2019

@author: zhangxuan
"""

import sys
import os
import datetime
import gffutils
import numpy as np

# 读gtf文件
db = gffutils.create_db("transcript.gtf", dbfn=":memory:")

print("Feature types = ", list(db.featuretypes()))

geneSubFeature = "transcript"

for feature_type in db.featuretypes():
    counts = db.count_features_of_type(featuretype=feature_type)
    print(feature_type, counts)

levelTwoType = []

geneDist = os.path.join(outdir, "gene.dist.tsv")
with open(geneDist, 'w') as f:
    f.write("id\ttranscript_num\tgc\tgc1\tgc2\tgc3\tlength\n")
    for gene in db.features_of_type("gene"):
        transcriptCounts = str(len(list(db.children(gene))))
        transcriptType = [t.featuretype for t in db.children(gene, level=1)]
        levelTwoType += transcriptType
        geneFa = gene.sequence(fasta)
        gc = GC(geneFa)
        gc123 = GC123(geneFa)
        geneLen = gene.end - gene.start + 1
        items = [gene.id, transcriptCounts, str(gc), str(gc123[1]), str(gc123[2]), str(gc123[3]), str(geneLen)]
        linestr = '\t'.join(items)
        f.write(linestr + '\n')


print(set(levelTwoType))

# 基因间区长度分布

chroms = [i['seqid'] for i in db.execute('SELECT DISTINCT seqid FROM features;')]

def genes_on_chrom(chroms):
    """
    Yield genes on `chrom`, sorted by start position
    """
    for g in db.features_of_type('gene', order_by='start', limit=(chroms, 0, 1e9)):
        g.strand = '.'
        yield g

def intergenic(chroms):
    """
    Yield intergenic features
    """
    for chrom in chroms:
        genes = genes_on_chrom(chrom)
        for intergenic in db.interfeatures(features = genes, new_featuretype="intergenic"):
            yield intergenic


intergenicDist = os.path.join(outdir, "intergenic.dist.tsv")
with open(intergenicDist, 'w') as f:
    f.write("id\tgc\tgc1\tgc2\tgc3\tlength\n")
    for intergenic in intergenic(chroms):
        if intergenic.end - intergenic.start < 10:
            continue
        intergenic_id = '-'.join(intergenic.attributes['ID'])
        sequence = intergenic.sequence(fasta)
        gc = GC(sequence)
        gc123 = GC123(sequence)
        intergenicLen = intergenic.end - intergenic.start + 1
        items = [intergenic_id, str(gc), str(gc123[1]), str(gc123[2]), str(gc123[3]), str(intergenicLen)]
        linestr = '\t'.join(items)
        f.write(linestr + '\n')

for type in list(set(levelTwoType)):
    rnaDist = os.path.join(outdir, "{}.dist.tsv".format(type))
    exonDist = os.path.join(outdir, "{}.exon.dist.tsv".format(type))
    intronDist = os.path.join(outdir, "{}.intron.dist.tsv".format(type))
    cdsDist = os.path.join(outdir, "{}.cds.dist.tsv".format(type))
    rna_f = open(rnaDist, 'w')
    exon_f = open(exonDist, 'w')
    intron_f = open(intronDist, 'w')
    cds_f = open(cdsDist, 'w')
    rna_f.write("id\texon_num\tgc\tgc1\tgc2\tgc3\tlength\n")
    exon_f.write("id\tgc\tgc1\tgc2\tgc3\tlength\n")
    intron_f.write("id\tgc\tgc1\tgc2\tgc3\tlength\n")
    cds_f.write("id\tgc\tgc1\tgc2\tgc3\tlength\n")
    for gene in db.features_of_type('gene'):
        for transcript in db.children(gene, featuretype=type):
            exonCounts = str(len(list(db.children(transcript))))
            sequence = transcript.sequence(fasta)
            if not sequence:
                    sys.stderr.write(transcript.id + 'failed\n')
                    continue
            gc = GC(sequence)
            gc123 = GC123(sequence)
            transcriptLen = transcript.end - transcript.start + 1
            items = [transcript.id, str(gc), str(gc123[1]), str(gc123[2]), str(gc123[3]), str(intergenicLen)]
            linestr = '\t'.join(items)
            rna_f.write(linestr + '\n')
            for exon in db.children(transcript, featuretype='exon'):
                sequence = exon.sequence(fasta)
                if not sequence:
                    sys.stderr.write(exon.id + 'failed\n')
                    continue
                gc = GC(sequence)
                gc123 = GC123(sequence)
                exonLen = exon.end - exon.start + 1
                items = [exon.id, exonCounts, str(gc), str(gc123[1]), str(gc123[2]), str(gc123[3]), str(exonLen)]
                linestr = '\t'.join(items)
                exon_f.write(linestr + '\n')
            for intron in db.children(transcript, featuretype='intron'):
                sequence = intron.sequence(fasta)
                if not sequence:
                    sys.stderr.write(intron.id + 'failed\n')
                    continue
                gc = GC(sequence)
                gc123 = GC123(sequence)
                intronLen = intron.end - intron.start + 1
                items = [intron.id, str(gc), str(gc123[1]), str(gc123[2]), str(gc123[3]), str(intronLen)]
                linestr = '\t'.join(items)
                intron_f.write(linestr + '\n')
            for cds in db.children(transcript, featuretype='CDS'):
                sequence = cds.sequence(fasta)
                if not sequence:
                    sys.stderr.write(cds.id + 'failed\n')
                    continue
                gc = GC(sequence)
                gc123 = GC123(sequence)
                cdsLen = cds.end - cds.start + 1
                items = [cds.id, str(gc), str(gc123[1]), str(gc123[2]), str(gc123[3]), str(cdsLen)]
                linestr = '\t'.join(items)
                cds_f.write(linestr + '\n')
    
    rna_f.close()
    exon_f.close()
    intron_f.close()
    cds_f.close()