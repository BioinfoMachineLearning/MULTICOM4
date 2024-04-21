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
        with open(pkl, 'rb') as f:
            prediction_result = pickle.load(f)
            print(prediction_result['ranking_confidence'])
            ranking_confidences[pdb.replace('unrelaxed_', '').replace('.pdb', '')] = float(prediction_result['ranking_confidence'])
            ranked_order.append(pdb.replace('unrelaxed_', '').replace('.pdb', ''))

    # Rank by model confidence and write out relaxed PDBs in rank order.
    ranked_order = []
    for idx, (model_name, _) in enumerate(
            sorted(ranking_confidences.items(), key=lambda x: x[1], reverse=True)):
        ranked_order.append(model_name)
        ranked_output_path = os.path.join(args.indir, f'ranked_{idx}.pdb')
        srcpdb = f"relaxed_{model_name}.pdb"
        if not os.path.exists(srcpdb):
            srcpdb = f"unrelaxed_{model_name}.pdb"
        os.system(f"cp {srcpdb} {ranked_output_path}")

    ranking_output_path = os.path.join(args.indir, 'ranking_debug.json')
    with open(ranking_output_path, 'w') as f:
        label = 'iptm+ptm' if 'iptm' in prediction_result else 'plddts'
        f.write(json.dumps(
            {label: ranking_confidences, 'order': ranked_order}, indent=4))
