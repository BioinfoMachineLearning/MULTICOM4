import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from multicom_dev.common.util import is_file, is_dir, makedir_if_not_exists, clean_dir

pdb2seq = '/home/bml_casp15/BML_CASP15/utils/pdb2seq_chain_all.pl'
PDB_CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=is_dir, required=True)
    parser.add_argument('--outdir', type=str, required=True)

    args = parser.parse_args()

    for target in os.listdir(args.indir):
        all_chain_seq = ''
        for model in os.listdir(args.indir + '/' + target):
            chain_seqs = os.popen(f"perl {pdb2seq} {args.indir}/{target}/{model}").readlines()
            if len(all_chain_seq) == 0:
                all_chain_seq = '_'.join(chain_seqs)
            elif all_chain_seq != '_'.join(chain_seqs):
                raise Exception("The sequence are not the same!")

        with open(args.outdir + '/' + target + '.fasta', 'w') as fw:
            seqs = all_chain_seq.split('_')
            for i in range(len(seqs)):
                fw.write(f">{target}_{PDB_CHAIN_IDS[i]}\n{seqs[i]}")

