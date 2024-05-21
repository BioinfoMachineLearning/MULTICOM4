import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle

def run_commands(inparams):
    cmd = inparams[0]
    os.system(cmd)

class Gate_qa:

    def __init__(self, params):
        self.params = params

    def run_monomer_qa(self, fasta_path, input_dir, outputdir, contact_map_file = "", dist_map_file = ""):

        resultfile = os.path.join(outputdir, 'ensemble_af.csv')

        if not os.path.exists(resultfile):

            process_list = []

            if not os.path.exists(contact_map_file):
                cmd = f"sh {self.params['dncon4_program']} {fasta_path} {outputdir}/dncon4 &> {outputdir}/dncon4.log"
                process_list.append([cmd])

            if not os.path.exists(dist_map_file):
                cmd = f"sh {self.params['deepdist_program']} {fasta_path} {outputdir}/deepdist &> {outputdir}/deepdist.log"
                process_list.append([cmd])

            if len(process_list) > 0:
                pool = Pool(processes=len(process_list))
                results = pool.map(run_commands, process_list)
                pool.close()
                pool.join()

            targetname = open(fasta_path).readlines()[0].rstrip('\n')[1:]
            if not os.path.exists(contact_map_file):
                contact_map_file = os.path.join(outputdir, 'dncon4', f'{targetname}.dncon2.rr')
            
            if not os.path.exists(dist_map_file):
                dist_map_file = os.path.join(outputdir, 'deepdist', f'{targetname}.txt')

            cmd = f"sh {self.params['gate_qa_program_dir']} monomer {fasta_path} {input_dir} {outputdir} {contact_map_file} {dist_map_file}"
            os.system(cmd)

            if not os.path.exists(resultfile):
                raise Exception(f"Failed to run gate qa!")

        return resultfile

    def run_multimer_qa(self, fasta_path, input_dir, pkl_dir, outputdir, use_af_features=True):
                            
        resultfile = os.path.join(outputdir, 'ensemble_af.csv')
        if not use_af_features:
            resultfile = os.path.join(outputdir, 'ensemble_nonaf.csv')

        if os.path.exists(resultfile) and len(pd.read_csv(resultfile)) == len(os.listdir(input_dir)):
            return resultfile

        cmd = f"sh {self.params['gate_qa_program_dir']} multimer {fasta_path} {input_dir} {outputdir} {pkl_dir} "
        cmd += "True" if use_af_features else "False"
        os.system(cmd)

        return resultfile
