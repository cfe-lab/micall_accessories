import argparse
from micall.core.denovo import genotype

def main(args):
    with open('./blast.csv', 'w') as o:
        genotype(fasta=args.fasta, blast_csv=o, group_refs={})

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta')

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)
