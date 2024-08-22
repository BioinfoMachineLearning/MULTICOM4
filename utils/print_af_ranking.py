import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
from multicom4.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file, is_file, is_dir
from multicom4.multimer_structure_evaluation.alphafold_ranking import Alphafold_pkl_qa
import json

PDB_CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
def _reorder_chains(pdbstring):
    new_pdbstring = []
    first_chain_id = None
    for line in pdbstring.split('\n'):
        if line.startswith('ATOM') or line.startswith('TER'):
            chain_id = line[21]
            if first_chain_id is None:
                first_chain_id = chain_id
            if first_chain_id != "A":
                new_pdbstring += [line[:21] + PDB_CHAIN_IDS[PDB_CHAIN_IDS.find(chain_id)-1] + line[22:]]
            else:
                new_pdbstring += [line]
        else:
            new_pdbstring += [line]
    return '\n'.join(new_pdbstring)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=is_file, required=True)
    args = parser.parse_args()

    ranking_confidences = {}
    ranked_order = []

    for pdb in os.listdir(args.indir):
        if pdb.find('unrelaxed') != 0:
            continue
        pkl = pdb.replace('unrelaxed', 'result').replace('.pdb', '.pkl')
        with open(args.indir + '/' + pkl, 'rb') as f:
            prediction_result = pickle.load(f)
            print(prediction_result['ranking_confidence'])
            ranking_confidences[pdb.replace('unrelaxed_', '').replace('.pdb', '')] = float(prediction_result['ranking_confidence'])
            ranked_order.append(pdb.replace('unrelaxed_', '').replace('.pdb', ''))

    # Rank by model confidence and write out relaxed PDBs in rank order.
    ranked_order = []
    for idx, (model_name, _) in enumerate(
            sorted(ranking_confidences.items(), key=lambda x: x[1], reverse=True)):
        print(f"{model_name}\t{ranking_confidences[model_name]}")

