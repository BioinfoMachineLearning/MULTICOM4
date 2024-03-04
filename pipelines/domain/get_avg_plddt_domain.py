import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import json
import numpy as np

def find_idx_in_domain_ranges(residue_idx, domain_starts, domain_ends):
    for start, end in zip(domain_starts, domain_ends):
        if int(start) <= residue_idx <= int(end):
            return True
    return False

def get_avg_factor(infile, domain_starts, domain_ends):
    bfactors = []
    read = set()
    for line in open(infile):
        if line.startswith('ATOM'):
            resn = line[17:20]
            chain = line[21]
            resi = int(line[22:26])

            if not find_idx_in_domain_ranges(resi, domain_starts, domain_ends):
                continue

            icode = line[26]
            bfactor = line[60:66]
            r_uid = (resn, chain, resi, icode)
            if r_uid not in read:
                read.add(r_uid)
                #print(read)
            else:
                continue
            bfactors += [float(bfactor)]
    return np.mean(np.array(bfactors))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--domain_info', type=str, required=True)
    parser.add_argument('--workdir', type=str, required=True)
    args = parser.parse_args()

    all_domain_ranges = []

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
            # list index start from 0
            starts += [int(start)]
            ends += [int(end)]

        all_domain_ranges += [(starts, ends)]

    alphafold_ranking_file = os.path.join(args.workdir, 'ranking_debug.json')
    ranking_json = json.loads(open(alphafold_ranking_file).read())

    for model_name in ranking_json["plddts"]:
        line = [model_name]
        predict_model = os.path.join(args.workdir, f"relaxed_{model_name}.pdb")
        for domain_idx, (domain_starts, domain_ends) in enumerate(all_domain_ranges):
            domain_plddt = get_avg_factor(predict_model, domain_starts, domain_ends)
            line += [str(domain_plddt)]
        print(' '.join(line))
