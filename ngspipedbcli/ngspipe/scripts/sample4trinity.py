import os
import sys

sampledir = sys.argv[1]

samples = os.listdir(sampledir)
samples = sorted(samples)

for sample in samples:
    expr = sample.split("-")
    print(expr[0], sample, os.path.join(sampledir,sample,sample+".1.non_rRNA_val_1.fq"), os.path.join(sampledir,sample,sample+".2.non_rRNA_val_2.fq"))
