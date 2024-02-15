import re, os, math
from collections import OrderedDict, defaultdict
from copy import deepcopy
import numpy as np
import pandas as pd
from multicom4.monomer_alignment_generation.alignment import *
from multicom4.monomer_alignments_concatenation.species_interact_v3 import Species_interact_v3
import itertools

def combineMSAs(MSA_Tags, is_homomers=False):  # dict

    combined_MSA_Tags = []

    if is_homomers:
        num_examples = len(list(MSA_Tags.keys()))
        for chain in MSA_Tags:
            for msa_name in MSA_Tags[chain]:
                combined_MSA_Tags += [[msa_name] * num_examples]
            break
    else:
        x = []
        for chain in MSA_Tags:
            x.append(MSA_Tags[chain])
        combined_MSA_Tags = list(itertools.product(*x))

    return combined_MSA_Tags


class DeepMSA2_pairing:

    def get_pairs(deepmsa_chain_alignments, is_homomers=False):
        num_examples = len(list(deepmsa_chain_alignments.keys()))
        MSA_Tags = {}
        for chain in deepmsa_chain_alignments:
            MSA_Tags[chain] = list(deepmsa_chain_alignments[chain].keys())

        combined_MSA_Tags = combineMSAs(MSA_Tags, is_homomers)

        combined_MSA_alignments = {}
        for combined_MSA_Tag in combined_MSA_Tags:  # c(qMSA,DeepMSA.hhb,DeepJGI3)
            MSA_name = '_'.join([combined_MSA_Tag[i] for i in range(num_examples)])
            deepmsa_alignments = []
            for i, chain in enumerate(deepmsa_chain_alignments):
                deepmsa_alignments += [deepmsa_chain_alignments[chain][combined_MSA_Tag[i]]]
            combined_MSA_alignments[MSA_name] = deepmsa_alignments

        return combined_MSA_alignments

    def rank_msas(deepmsa_complex_aln_dir, deepmsa_ranking_files, outfile, calNf):

        deepmsa_ranking_info = []
        for deepmsa_ranking_file in deepmsa_ranking_files:

            deepmsa_ranking = {}
            contents = open(deepmsa_ranking_file).readlines()     

            for index, line in enumerate(contents):
                line = line.rstrip('\n')
                msa_name, plddt = line.split()
                msa_name = msa_name.replace('_', '.')
                deepmsa_ranking[msa_name] = float(plddt)
            deepmsa_ranking_info += [deepmsa_ranking]

        data_dict = {'name': [], 'fscore': [], 'neff': []}
        for combination in os.listdir(deepmsa_complex_aln_dir):

            data_dict['name'] += [combination]
            print(f"Calculate Neff for joint MSA {combination}")

            joint_MSA_path = os.path.join(deepmsa_complex_aln_dir, combination, combination + '.a3m')

            joint_aln_path = os.path.join(deepmsa_complex_aln_dir, combination, combination + '.aln')
            contents = open(joint_MSA_path).readlines()
            with open(joint_aln_path, 'w') as fw:
                for line in contents:
                    if line[0] == ">":
                        continue
                    fw.write(line)

            Neff_MSA_path= joint_MSA_path + ".neff"

            cmd = f"{calNf} {joint_aln_path} > {Neff_MSA_path}"
            print(cmd)
            os.system(cmd)

            #### read neff 
            contents = open(Neff_MSA_path).readlines()

            neff = float(contents[0].rstrip('\n'))

            msa_names = combination.split('_')
            Nchain = len(msa_names)
            fscore = np.mean(np.array([deepmsa_ranking_info[i][msa_names[i]] for i in range(Nchain)]))
            
            fscore = math.log(neff, 10) * fscore

            data_dict['fscore'] += [fscore]
            data_dict['neff'] += [neff]

        df = pd.DataFrame(data_dict)
        df = df.sort_values(by='fscore', ascending=False)
        df.reset_index(inplace=True, drop=True)
        df.to_csv(outfile)
        

