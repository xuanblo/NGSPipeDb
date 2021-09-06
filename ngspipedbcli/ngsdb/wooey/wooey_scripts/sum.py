import argparse
import sys

parser = argparse.ArgumentParser(description="Find the sum of all the numbers below a certain number.")
parser.add_argument('--input', help='The number to find the sum of numbers below.', type=str, default='mm.fa')
parser.add_argument('--outfile', help='The number to find the sum of numbers below.', type=str, default='mm.fa')

def main():
    args = parser.parse_args()
    i = args.input
    with open(args.outfile, 'w') as f:
        f.write(i)
    return 0

if __name__ == "__main__":
    sys.exit(main())
