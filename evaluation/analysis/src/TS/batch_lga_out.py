import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd

def run_command(inparams):
    lga_program, pdb = inparams
    cmd = f"{lga_program} -3 -sda -ch2:A {pdb}"
    print(cmd)
    os.system(cmd)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str, default='MOL2')
    parser.add_argument('--outdir', type=str, default='TMP')
    parser.add_argument('--lga_program', type=str, required=True)

    args = parser.parse_args()

    process_list = []

    for targetname in os.listdir(args.indir):
        
        outdir = args.outdir + '/' + targetname

        os.makedirs(outdir, exist_ok=True)

        for model in os.listdir(args.indir + '/' + targetname):
            
            result_file = os.path.join(args.outdir, targetname + '/' + model + '.lga')
            if os.path.exists(result_file) and len(open(result_file).readlines()) >= 5:
                continue

            process_list.append([args.lga_program, targetname + '/' + model])

    pool = Pool(processes=120)
    results = pool.map(run_command, process_list)
    pool.close()
    pool.join()
