import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_path', type=str, required=True)
    parser.add_argument('--domain_info', type=str, required=True)
    parser.add_argument('--outputdir', type=str, required=True)
    args = parser.parse_args()

    os.makedirs(args.outputdir, exist_ok=True)

    sequence = open(args.fasta_path).readlines()[1].rstrip('\n')
    for line in open(args.domain_info):
        line = line.rstrip('\n')
        # domain 1: 140-317 Normal
        # domain 1: 140-317
        domain_name, domain_range = line.split(':')
        domain_name = domain_name.replace(' ', '')
        domain_range = domain_range.split()[0]

        domain_sequence = ""
        starts, ends = [], []

        for domain_range_str in domain_range.split(','):
            start, end = domain_range_str.split('-')
            starts += [int(start)-1]
            ends += [int(end)]

        for start, end in zip(starts, ends):
            domain_sequence += sequence[start:end]

        with open(os.path.join(args.outputdir, domain_name + '.fasta'), 'w') as fw:
            fw.write(f">{domain_name}\n{domain_sequence}")
