import re
from collections import OrderedDict, defaultdict
from copy import deepcopy
import numpy as np
import pandas as pd
from bml_casp15.monomer_alignment_generation.alignment import *
import os

MSA_CROP_SIZE = 1024

def read_fasta(infasta):
    targets = []
    seqs = []
    for line in open(infasta):
        line = line.rstrip('\n')
        if line[0] == ">":
            targets += [line[1:]]
        elif len(line) > 0:
            seqs += [line]
    return targets, seqs


class Combine_alignments_v2:

    def __init__(self, params):
        self.params = params

    def combine_interactions(self, methods, outdir):
        
        pairs_alignments = {}

        chain_alignments = {}

        for method in methods:

            interact_csv = f"{outdir}/{method}/{method}_interact.csv"
            
            if not os.path.exists(interact_csv):
                raise Exception(f"Cannot find the interact csv {interact_csv} for {method}!")

            pairs = pd.read_csv(interact_csv)
            print(pairs)
            num_chains = int(pairs.shape[1] / 2)

            chain_seqs = {}
            chain_headers = {}
            print("11111111111111111111")
            for chain_index in range(num_chains):
                chain_name = pairs.loc[0, f"id_{chain_index+1}"]
                names, seqs = read_fasta(f"{outdir}/{method}/{chain_name}_con.a3m")
                chain_seqs[chain_index] = seqs
                chain_headers[chain_index] = names
                if chain_index not in chain_alignments:
                    chain_alignments[chain_index] = dict(names=[], seqs=[], seq_len=len(seqs[0]), chain_name=chain_name, chain_seq=seqs[0])

            print("222222222222222222222")
            for index in range(1, len(pairs)):
                pair_info = []
                pair_seqs = []
                pair_headers = []

                print("444444444444444444444444")
                for chain_index in range(num_chains):
                    if pairs.loc[index, f"id_{chain_index+1}"].find('placeholder') < 0:
                        pair_info += [str(chain_index)]
                        pair_seqs += [chain_seqs[chain_index][index]]
                        pair_headers += [chain_headers[chain_index][index]]
                
                print("55555555555555555555")
                if "_".join(pair_info) not in pairs_alignments:
                    pairs_alignments["_".join(pair_info)] = dict(seqs=[], headers=[])
                
                print("6666666666666666666666666666")
                pairs_alignments["_".join(pair_info)]['seqs'] += ["".join(pair_seqs)]
                pairs_alignments["_".join(pair_info)]['headers'] += ["_____".join(pair_headers)]

        outputdir = outdir + '/combine_v2'
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)
        
        # do hhfilter
        total_msa_depth = np.sum(np.array([len(pairs_alignments[pair_info]['seqs']) for pair_info in pairs_alignments]))
        pairs_alignments_filt = {}
        if total_msa_depth <= MSA_CROP_SIZE:
            pairs_alignments_filt = deepcopy(pairs_alignments)

        else:
            total_pairs = len([1 for pair_info in pairs_alignments])
        
            per_pair_count = int(MSA_CROP_SIZE / total_pairs) 
            for pair_info in pairs_alignments:
                outfile = outputdir + '/' + pair_info + '.a3m'
                with open(outfile, 'w') as fw:
                    for i in range(len(pairs_alignments[pair_info]['seqs'])):
                        name = pairs_alignments[pair_info]['headers'][i]
                        seq = pairs_alignments[pair_info]['seqs'][i]
                        fw.write(f">{name}\n{seq}\n")

                if len(pairs_alignments[pair_info]['seqs']) > per_pair_count:
                    cmd = f"{self.params['hhfilter_program']} -diff {per_pair_count} -i {outfile} -o {outfile}.filt -id 80 -cov 50"
                    os.system(cmd)
                    names, seqs = read_fasta(outfile + '.filt')
                else:
                    names, seqs = read_fasta(outfile)

                pairs_alignments_filt[pair_info] = dict(seqs=[], headers=[])
                for name, seq in zip(names, seqs):
                    pairs_alignments_filt[pair_info]['headers'] += [name]
                    pairs_alignments_filt[pair_info]['seqs'] += [seq]
        
        print("33333333333333333333333")

        # recombine all the alignments
        count = 0
        for pair_info in pairs_alignments_filt:
            names = pairs_alignments_filt[pair_info]['headers']
            seqs = pairs_alignments_filt[pair_info]['seqs']

            for name, seq in zip(names, seqs):
                headers = name.split('_____')
                seq_start_idx = 0
                for i in chain_alignments:
                    if str(i) not in pair_info.split('_'):
                        chain_alignments[i]['names'] += ['placeholder' + str(count)]
                        chain_alignments[i]['seqs'] += ['-' * int(chain_alignments[i]['seq_len'])]
                    else:
                        idx_count = pair_info.split('_').index(str(i))
                        chain_alignments[i]['names'] += [headers[idx_count]]
                        chain_alignments[i]['seqs'] += [seq[seq_start_idx:seq_start_idx+int(chain_alignments[i]['seq_len'])]]
                        seq_start_idx += int(chain_alignments[i]['seq_len'])
                count += 1

        fw = open(outputdir + '/combine.a3m', 'w')
        depth = 0
        chain_names = []
        chain_seqs = []
        filter_pair_ids = {}
        for chain_idx in chain_alignments:
            if depth == 0:
                depth = len(chain_alignments[chain_idx]['seqs'])
            elif depth != len(chain_alignments[chain_idx]['seqs']):
                raise Exception("The msa depth for each chains are different!")
            chain_alignments[chain_idx]['fw'] = open(outputdir + '/' + chain_alignments[chain_idx]['chain_name'] + '_con.a3m', 'w')
            chain_alignments[chain_idx]['fw'].write(f">{chain_alignments[chain_idx]['chain_name']}\n{chain_alignments[chain_idx]['chain_seq']}\n")
            chain_seqs += [chain_alignments[chain_idx]['chain_seq']]
            chain_names += [chain_alignments[chain_idx]['chain_name']]
            filter_pair_ids[f"id_{chain_idx + 1}"] = [chain_alignments[chain_idx]['chain_name']]
            filter_pair_ids[f"index_{chain_idx + 1}"] = [0]

        fw.write(f">{'_____'.join(chain_names)}\n{''.join(chain_seqs)}\n")

        
        for i in range(depth):
            headers = []
            seqs = []
            for chain_idx in chain_alignments:
                chain_alignments[chain_idx]['fw'].write(f">{chain_alignments[chain_idx]['names'][i]}\n{chain_alignments[chain_idx]['seqs'][i]}\n")
                headers += [chain_alignments[chain_idx]['names'][i]]
                seqs += [chain_alignments[chain_idx]['seqs'][i]]

                filter_pair_ids[f"id_{chain_idx + 1}"] += [chain_alignments[chain_idx]['names'][i]]
                filter_pair_ids[f"index_{chain_idx + 1}"] += [i+1]

            fw.write(f">{'_____'.join(headers)}\n{''.join(seqs)}\n")
        
        fw.close()
        for chain_idx in chain_alignments:
            chain_alignments[chain_idx]['fw'].close()

        pd.DataFrame(filter_pair_ids).to_csv(f"{outputdir}/combine_v2_interact.csv", index=False)