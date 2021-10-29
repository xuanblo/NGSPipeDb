import os
import sys
import argparse
from Bio import SeqIO
import pprint

parser = argparse.ArgumentParser(description='get sub genome from given parameters')
parser.add_argument('-g', '--genome', required=True, help='genome fasta file')
parser.add_argument('-go', '--subgenomeout', required=True, help='sub genome fasta file')
parser.add_argument('-r', '--region', required=True, help='all (whole genome) or chr1:10000-20000 (seq name:start position-end position) or chr1 (seq name)')

args = parser.parse_args()

# fasta
if args.region == 'all':
    sub_fasta_command = 'cp {} {}'.format(args.genome, args.subgenomeout)
    os.system(sub_fasta_command)
else:
    record_dict = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
    width = 60
    gout = open(args.subgenomeout, 'w')
    if ':' in args.region: # interval
        fid, pos = args.region.split(':')
        start ,end = [int(i)-1 for i in pos.split('-')]
        gout.write('>'+args.region+'\n')  # use any record ID
        seq_str = str(record_dict[fid].seq[start:end])
        new_seq_str = [seq_str[i:i+width] for i in range(0, len(seq_str), width)]
        gout.write('\n'.join(new_seq_str)+'\n')  # use any record ID
    else: # id
        gout.write('>'+record_dict[args.region].id+'\n')  # use any record ID
        seq_str = str(record_dict[args.region].seq)
        new_seq_str = [seq_str[i:i+width] for i in range(0, len(seq_str), width)]
        gout.write('\n'.join(new_seq_str)+'\n')  # use any record ID

