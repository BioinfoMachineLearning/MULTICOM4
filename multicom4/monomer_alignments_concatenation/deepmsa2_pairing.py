import re, os, math
from collections import OrderedDict, defaultdict
from copy import deepcopy
import numpy as np
import pandas as pd
from multicom4.monomer_alignment_generation.alignment import *
from multicom4.monomer_alignments_concatenation.species_interact_v3 import Species_interact_v3
import itertools
from multicom4.common import config

def getTopN(maxP, chainN, is_homomers):
    topN = 10
    if is_homomers:
        topN = 50
    else:
        topN = math.floor(math.pow(maxP, 1.0 / chainN))
    return topN

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


class DeepMSA2_pairing(config.pipeline):

    def __init__(self, is_homomers=False):

        super().__init__()

        self.deepmsa2_config = self.homomer_config.predictors.deepmsa2 if is_homomers else self.heteromer_config.predictors.deepmsa2

        self.is_homomers = is_homomers

    def get_pairs(self, chain_alignments):

        unique_seuqences = {}

        for chain in chain_alignments:
            if chain == "outdir":
                continue
            chain_sequence = chain_alignments[chain]['chain_seq']
            if chain_sequence not in unique_seuqences:
                unique_seuqences[chain_sequence] = []
            unique_seuqences[chain_sequence] += [chain]

        topN = getTopN(self.deepmsa2_config.max_pairs, len(unique_seuqences), self.is_homomers)

        MSA_Tags = {}
        for sequence in unique_seuqences:
            MSA_Tags[sequence] = []
            deepmsa_ranking_file = chain_alignments[unique_seuqences[sequence][0]]['deepmsa_ranking_file']
            contents = open(deepmsa_ranking_file).readlines()
            for index, line in enumerate(contents):
                if index >= topN:
                    break
                line = line.rstrip('\n')
                msa_name, plddt = line.split()
                msa_name = msa_name.replace('_', '.')
                MSA_Tags[sequence] += [msa_name]

        combined_MSA_Tags = combineMSAs(MSA_Tags, self.is_homomers)

        combined_MSA_alignments = {}
        for combined_MSA_Tag in combined_MSA_Tags:  # c(qMSA,DeepMSA.hhb,DeepJGI3)
            MSA_names = []
            deepmsa_alignments = []
            for i, sequence in enumerate(unique_seuqences):
                aln_name = combined_MSA_Tag[i]
                for chain_id in unique_seuqences[sequence]:
                    with open(chain_alignments[chain_id][aln_name]) as f:
                        deepmsa_alignment = Alignment.from_file(f, format="a3m", a3m_inserts="delete")
                    deepmsa_alignments += [deepmsa_alignment]
                    MSA_names += [aln_name]

            combined_MSA_alignments['_'.join(MSA_names)] = deepmsa_alignments

        return combined_MSA_alignments

    def rank_msas(self, deepmsa_complex_aln_dir, deepmsa_ranking_files, outfile, calNf):

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
        

