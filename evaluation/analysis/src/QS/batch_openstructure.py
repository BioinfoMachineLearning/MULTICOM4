import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from multicom4.common.util import is_file, is_dir, makedir_if_not_exists, clean_dir
import json

ema_exec = 'docker run --rm -v $(pwd):/home -u $(id -u ${USER}):$(id -g ${USER}) registry.scicore.unibas.ch/schwede/openstructure'

def run_command(inparams):
    mdl_path, trg_path, result_path = inparams

    cmd = [ema_exec, "compare-structures",
                       "-m", mdl_path, "-r", trg_path,
                       "-mf", "pdb", "-rf", "pdb", "-rna", "--out", result_path,
                       "--ics",
                       "--ips",
                       "--qs-score",
                       "--lddt",
                       "--local-lddt",
                       #"--cad-score",
                       #"--local-cad",
                       #"--rigid-scores", 
                       "--patch-scores",
                       "--tm-score",
                       "--dockq",
                       "--ics-trimmed",
                       "--ips-trimmed"]

    cmd = ' '.join(cmd)
    print(cmd)
    os.system(cmd)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str, required=True)
    parser.add_argument('--nativedir', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)

    args = parser.parse_args()

    process_list = []

    makedir_if_not_exists(args.outdir)

    for native_pdb in os.listdir(args.nativedir):

        targetname = native_pdb.rstrip('_filtered.pdb')
        
        if not os.path.exists(args.indir + '/' + targetname):
            continue

        outdir = args.outdir + '/' + targetname

        makedir_if_not_exists(outdir)

        for model in os.listdir(args.indir + '/' + targetname):
            
            outfile =  outdir + '/' + model + '_filtered_out'
            if os.path.exists(outfile):
                try:
                    with open(outfile) as f:
                        data = json.load(f)
                        if data['qs_best'] is not None and data['tm_score'] is not None:
                            continue
                except Exception as e:
                    
                    print(e)

                    os.system(f"cp {args.indir}/{targetname}/{model} {outdir}/{model}_filtered")

                    #process_list.append([args.usalign_program, args.nativedir + '/' + native_pdb, args.indir + '/' + targetname + '/' + model, outdir + '/' + model + '_filtered_out'])
                    process_list.append([args.indir + '/' + targetname + '/' + model, args.nativedir + '/' + native_pdb, outdir + '/' + model + '_filtered_out'])
            
            else:

                os.system(f"cp {args.indir}/{targetname}/{model} {outdir}/{model}_filtered")

                #process_list.append([args.usalign_program, args.nativedir + '/' + native_pdb, args.indir + '/' + targetname + '/' + model, outdir + '/' + model + '_filtered_out'])
                process_list.append([args.indir + '/' + targetname + '/' + model, args.nativedir + '/' + native_pdb, outdir + '/' + model + '_filtered_out'])
            
                

    pool = Pool(processes=100)
    results = pool.map(run_command, process_list)
    pool.close()
    pool.join()
