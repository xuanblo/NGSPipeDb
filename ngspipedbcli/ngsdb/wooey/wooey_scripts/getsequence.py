import argparse
import sys

parser = argparse.ArgumentParser(description="get sequence from genome")
parser.add_argument('--pos', help='chr:start-end', type=str, default='chr1:10-100')
parser.add_argument('--strand', help='+/-', type=str, default='+')

def main():
    args = parser.parse_args()
    print('developling')
    return 0

if __name__ == "__main__":
    sys.exit(main())
