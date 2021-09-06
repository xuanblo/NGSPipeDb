import argparse
import sys

parser = argparse.ArgumentParser(description="convert dna sequences to rna sequences.")
parser.add_argument('--input', help='dna sequences.', type=str, default='ATGC')

def main():
    args = parser.parse_args()
    i = args.input
    i = i.upper()
    print(i.replace('T', 'U'))
    return 0

if __name__ == "__main__":
    sys.exit(main())
