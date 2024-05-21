import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle, json

def get_bfactors(infile):
    bfactors = []
    read = set()
    for line in open(infile):
        if line.startswith('ATOM'):
            resn = line[17:20]
            chain = line[21]
            resi = line[22:26]
            icode = line[26]
            bfactor = line[60:66]
            r_uid = (resn, chain, resi, icode)
            if r_uid not in read:
                read.add(r_uid)
                #print(read)
            else:
                continue
            bfactors += [float(bfactor)]
    return np.array(bfactors)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str, required=True)
    args = parser.parse_args()

    ranking_confidences = {}
    ranked_order = []
    for pdb in ['4', '10', '50', '80', '100']:
        pdbfile = pdb + '.pdb'
        plddts = get_bfactors(os.path.join(args.indir, pdbfile))

        result_output_path = os.path.join(args.indir, f'result_{pdb}.pkl')
        with open(result_output_path, 'wb') as f:
            pickle.dump({'plddt': plddts}, f, protocol=4)

        ranked_order.append(pdb)
        ranking_confidences[pdb] = np.mean(plddts)

    ranked_order = [model_name for model_name, confidence in sorted(ranking_confidences.items(), key=lambda x: x[1], reverse=True)]

    # Write out relaxed PDBs in rank order.
    for idx, model_name in enumerate(ranked_order):
        ranked_output_path = os.path.join(args.indir, f'ranked_{idx}.pdb')
        os.system(f"cp {os.path.join(args.indir, model_name + '.pdb')} {ranked_output_path}")

    ranking_output_path = os.path.join(args.indir, 'ranking_debug.json')
    with open(ranking_output_path, 'w') as f:
        f.write(json.dumps({'plddts': ranking_confidences, 'order': ranked_order}, indent=4))


