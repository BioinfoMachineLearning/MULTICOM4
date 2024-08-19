import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
from multicom4.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file, is_file, is_dir

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    parser.add_argument('--mapping', type=str, required=True)
    args = parser.parse_args()

    for inpdb in os.listdir(args.indir):
        print(f"Processing {inpdb}")
        chain_contents = {}
        for line in open(args.indir + '/' + inpdb):
            if line.startswith('ATOM'):
                chain_id = line[21]
                if chain_id not in chain_contents:
                    chain_contents[chain_id] = [line]
                else:
                    chain_contents[chain_id] += [line]

        mappings = args.mapping.split(',')
        map_chains = {mapping.split('_')[0]:mapping.split('_')[1] for mapping in mappings}
        with open(args.outdir + '/' + inpdb, 'w') as fw:
            for chain_id in chain_contents:
                for line in chain_contents[chain_id]:
                    chain_trg = chain_id
                    if chain_id in map_chains:
                        chain_trg = map_chains[chain_id]
                    fw.write(line[:21] + chain_trg + line[22:])
                fw.write("TER\n")
