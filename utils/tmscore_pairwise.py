import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from multicom_dev.common.util import is_file, is_dir, makedir_if_not_exists, clean_dir
from multicom_dev.quaternary_structure_evaluation.pairwise_mmalign import *

def run_command(inparams):
    tmscore_program, input_dir, pdb1, pdb2, workdir = inparams
    cmd = f"{tmscore_program} {input_dir}/{pdb1} {input_dir}/{pdb2} > {workdir}/{pdb1}_{pdb2}"
    os.system(cmd)

    for line in open(f"{workdir}/{pdb1}_{pdb2}"):
        line = line.rstrip('\n')
        if len(line) == 0:
            continue
        contents = line.split()
        if contents[0] == 'TM-score':
            tmscore = float(contents[2])
        if contents[0] == 'GDT-score':
            gdtscore = float(contents[2])
            
    return pdb1, pdb2, tmscore
    
def run_pairwise(tmscore_program, indir, workdir, outfile):

    if os.path.exists(outfile):
        return

    ranking_pd = pd.DataFrame(columns=['model', 'avg_tmscore', 'avg_gdtscore'])

    pdbs = os.listdir(indir)

    process_list = []
    for i in range(len(pdbs)):
        for j in range(len(pdbs)):
            pdb1 = pdbs[i]
            pdb2 = pdbs[j]
            if pdb1 == pdb2:
                continue
            process_list.append([tmscore_program, indir, pdb1, pdb2, workdir])

    pool = Pool(processes=40)
    results = pool.map(run_command, process_list)
    pool.close()
    pool.join()

    tmscores_dict = {}
    gdtscores_dict = {}
    for result in results:
        pdb1, pdb2, tmscore, gdtscore = result
        tmscores_dict[f"{pdb1}_{pdb2}"] = tmscore
        gdtscores_dict[f"{pdb1}_{pdb2}"] = gdtscore

    for i in range(len(pdbs)):
        pdb1 = pdbs[i]
        ranking = {'model': pdb1}
        tmscores, gdtscores = [], []
        for pdb2 in pdbs:
            if pdb1 == pdb2:
                continue
            tmscores += [tmscores_dict[f"{pdb1}_{pdb2}"]]
            gdtscores += [gdtscores_dict[f"{pdb1}_{pdb2}"]]

        ranking['avg_tmscore'] = np.mean(np.array(tmscores))
        ranking['avg_gdtscore'] = np.mean(np.array(gdtscores))

        ranking_pd = ranking_pd.append(pd.DataFrame(ranking, index=[i]))

    ranking_pd.sort_values(by=['avg_tmscore'], ascending=False, ignore_index=True).to_csv(outfile)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=is_dir, required=True)
    parser.add_argument('--tmscore', type=is_file, required=True)
    parser.add_argument('--outdir', type=is_file, required=True)

    args = parser.parse_args()

    tempdir = args.outdir + '/temp'
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)

    for target in os.listdir(args.indir):
        workdir = tempdir + '/' + target
        if not os.path.exists(workdir):
            os.makedirs(workdir)
        run_pairwise(args.tmscore, args.indir + '/' + target, workdir, args.outdir + '/' + target + '.csv')
