import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd


def run_command(inparams):
    mmalign_program, input_dir, pdb1, pdb2, scorefile = inparams
    # cmd = mmalign_program + ' ' + os.path.join(input_dir, pdb1) + ' ' + os.path.join(input_dir, pdb2) + " | grep TM-score | awk '{print $2}' "
    # print(cmd)
    # tmscore_contents = os.popen(cmd).read().split('\n')
    # tmscore = float(tmscore_contents[1].rstrip('\n'))
    cmd = f"{mmalign_program} {input_dir}/{pdb1} {input_dir}/{pdb2} > {scorefile}"
    os.system(cmd)


def read_mmalign(infile):
    for line in open(infile):
        line = line.rstrip('\n')
        if len(line) > 0:
            if line.split()[0] == 'TM-score=' and line.find('Structure_2') > 0:
                tmscore = float(line.split()[1])
                return tmscore
    return 0


class Pairwise_MMalign_qa:

    def __init__(self, mmalign_program):

        self.mmalign_program = mmalign_program

    def run(self, input_dir, workdir):

        ranking_pd = pd.DataFrame(columns=['Name', 'MMalign score'])

        pdbs = sorted(os.listdir(input_dir))

        process_list = []
        for i in range(len(pdbs)):
            for j in range(len(pdbs)):
                pdb1 = pdbs[i]
                pdb2 = pdbs[j]
                if pdb1 == pdb2:
                    continue
                scorefile = os.path.join(workdir, f"{pdb1}_{pdb2}.mmalign")
                if not os.path.exists(scorefile) or len(open(scorefile).readlines()) < 5:
                    process_list.append([self.mmalign_program, input_dir, pdb1, pdb2, scorefile])

        pool = Pool(processes=60)
        results = pool.map(run_command, process_list)
        pool.close()
        pool.join()

        scores_dict = {}
        for i in range(len(pdbs)):
            for j in range(len(pdbs)):
                pdb1 = pdbs[i]
                pdb2 = pdbs[j]
                if pdb1 == pdb2:
                    continue
                scorefile = os.path.join(workdir, f"{pdb1}_{pdb2}.mmalign")
                if not os.path.exists(scorefile):
                    raise Exception(f"cannot find {scorefile}")
                scores_dict[f"{pdb1}_{pdb2}"] = read_mmalign(scorefile)

        for i in range(len(pdbs)):
            pdb1 = pdbs[i]
            ranking = {'Name': pdb1}
            scores = []
            for pdb2 in pdbs:
                if pdb1 == pdb2:
                    continue
                if f"{pdb1}_{pdb2}" in scores_dict:
                    score = scores_dict[f"{pdb1}_{pdb2}"]
                scores += [score]

            ranking['MMalign score'] = np.mean(np.array(scores))

            ranking_pd = ranking_pd.append(pd.DataFrame(ranking, index=[i]))

        return ranking_pd.sort_values(by=['MMalign score'], ascending=False, ignore_index=True)
