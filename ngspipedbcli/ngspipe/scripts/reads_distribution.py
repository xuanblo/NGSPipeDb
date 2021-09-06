import sys
import os
import pandas as pd

files = sys.argv[1:-1]

resultfile = sys.argv[-1]

samples = []
protein_coding = []
Others = []

for htseq_count in files:
    samplefile = os.path.basename(htseq_count)
    sample = samplefile.split('.')[0]
    protein_read_counts = 0
    other_read_counts = 0
    sys.stderr.write(htseq_count + '\n')
    with open(htseq_count, 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('__'):
                items = line.split('\t')
                protein_read_counts += int(items[1])
            else:
                if line.startswith('__no_feature'):
                    items = line.split('\t')
                    other_read_counts += int(items[1])
    samples.append(sample)
    total_read_counts = protein_read_counts + other_read_counts
    protein_coding.append("{}\n({}%)".format(protein_read_counts, round(protein_read_counts/total_read_counts, 2)))
    Others.append("{}\n({}%)".format(other_read_counts, round(protein_read_counts/total_read_counts, 2)))

data = {
    "Sample name": samples,
    "protein_coding": protein_coding,
    "Others": Others
}

df = pd.DataFrame(data)

df.to_csv(resultfile)