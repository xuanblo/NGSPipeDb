import argparse
import sys

parser = argparse.ArgumentParser(description="go enrichment analysis by clusterprofil")
parser.add_argument('--input', help='a list of differential expression genes', type=str, default='gene1,gene2')

def main():
    args = parser.parse_args()
    print('developling')
    return 0

if __name__ == "__main__":
    sys.exit(main())
