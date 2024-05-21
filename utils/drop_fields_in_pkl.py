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

    outdir = os.path.join(args.indir, 'reduced_pkls')
    os.makedirs(outdir ,exist_ok=True)

    for pkl in os.listdir(args.indir):
        if pkl.find('.pkl') < 0:
            continue
        if pkl == "features.pkl":
            continue

        result_output_path = os.path.join(outdir, pkl)

        with open(args.indir + '/' + pkl, 'rb') as f:
            prediction_result = pickle.load(f)

            keys_to_remove=['distogram', 'experimentally_resolved', 'masked_msa', 'aligned_confidence_probs']    
            d={}
            for k in prediction_result.keys():
                if k not in keys_to_remove:
                    d[k]=prediction_result[k]

            with open(result_output_path, 'wb') as f:
                pickle.dump(d, f, protocol=4)
            