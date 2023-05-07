from calendar import c
import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import random
import numpy as np
from Bio.PDB.PDBParser import PDBParser
from hungarian_algorithm import algorithm
import pandas as pd
from scipy.optimize import linear_sum_assignment

PDB_CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
dockq_program = '/home/bml_casp15/BML_CASP15/tools/DockQ/DockQ.py'
def parse_fasta(infasta):
    sequences = []
    descriptions = []
    index = -1
    for line in open(infasta):
        line = line.strip()
        if line.startswith('>'):
            index += 1
            descriptions.append(line[1:])  # Remove the '>' at the beginning.
            sequences.append('')
            continue
        elif not line:
            continue  # Skip blank lines.
        sequences[index] += line

    return sequences, descriptions

def make_chain_id_map(sequences, descriptions):
    """Makes a mapping from PDB-format chain ID to sequence and description."""
    if len(sequences) != len(descriptions):
        raise ValueError('sequences and descriptions must have equal length. '
                         f'Got {len(sequences)} != {len(descriptions)}.')
    chain_id_map = {}
    for chain_id, sequence, description in zip(
            PDB_CHAIN_IDS, sequences, descriptions):
        chain_id_map[chain_id] = sequence
    return chain_id_map


def split_pdb(complex_pdb, outdir):
    chain_pdbs = {}
    pre_chain = None
    i = 0
    for line in open(complex_pdb, 'r').readlines():
        if not line.startswith('ATOM'):
            continue
        chain_name = line[21]
        if pre_chain is None:
            pre_chain = chain_name
            fw = open(outdir + '/' + chain_name + '.pdb', 'w')
            chain_pdbs[chain_name] = {'pdb_path': outdir + '/' + chain_name + '.pdb', 'seq': ""}
            fw.write(line[:21] + ' ' + line[22:])
        elif chain_name == pre_chain:
            fw.write(line[:21] + ' ' + line[22:])
        else:
            fw.close()
            i = i + 1
            fw = open(outdir + '/' + chain_name + '.pdb', 'w')
            fw.write(line[:21] + ' ' + line[22:])
            pre_chain = chain_name
            chain_pdbs[chain_name] = {'pdb_path': outdir + '/' + chain_name + '.pdb', 'seq': ""}
    fw.close()

    for chain_id in chain_pdbs:
        chain_pdbs[chain_id]['seq'] = get_sequence(chain_pdbs[chain_id]['pdb_path'])

    return chain_pdbs


def cal_min_distance_between_chains(chain1, chain2, chain1_pdb, chain2_pdb):
    print(f"processing {chain1_pdb} and {chain2_pdb}")
    parser = PDBParser(PERMISSIVE=1)
    
    structure1 = parser.get_structure(chain1, chain1_pdb)
    structure2 = parser.get_structure(chain2, chain2_pdb)

    model1 = structure1[0]
    chain_id1 = list(model1.child_dict.keys())
    xyzPDB1 = model1[chain_id1[0]]

    model2 = structure2[0]
    chain_id2 = list(model2.child_dict.keys())
    xyzPDB2 = model2[chain_id2[0]]

    h_dist_map = np.zeros((len(xyzPDB1), len(xyzPDB2)))

    for i in range(1, len(xyzPDB1)+1):
        if i not in xyzPDB1:
            continue
        for j in range(1, len(xyzPDB2)+1):
            if j not in xyzPDB2:
                continue
            res_i = xyzPDB1[i]
            res_j = xyzPDB2[j]
            dist_list = []
            for atom_i in res_i:
                for atom_j in res_j:
                    if ('C' in atom_i.name or 'N' in atom_i.name or 'O' in atom_i.name or 'S' in atom_i.name) and \
                        ('C' in atom_j.name or 'N' in atom_j.name or 'O' in atom_j.name or 'S' in atom_j.name):
                        dist_list.append(atom_i - atom_j)
                    else:
                        continue
                    min_dist = np.min(dist_list)
                    h_dist_map[i-1, j-1] = min_dist

    min_dist = np.min(h_dist_map)
    return chain1, chain2, chain1_pdb, chain2_pdb, min_dist

def cal_dockq_score(inparams):
    key, inpdb, nativepdb = inparams
    cmd = f"python {dockq_program} -short " + inpdb + " " \
          + nativepdb + " | grep DockQ | awk '{print $2}'"
    score_default = 0.0
    try:
        score = os.popen(cmd).read()
        score = score.rstrip('\n')
        if len(score) > 0:
            score_default = float(score)
    except Exception as e:
        print(cmd)
        print(e)
    return key, score_default


def get_sequence(inpdb):
    """Enclosing logic in a function to simplify code"""

    seq_to_res_mapping = []
    res_codes = [
        # 20 canonical amino acids
        ('CYS', 'C'), ('ASP', 'D'), ('SER', 'S'), ('GLN', 'Q'),
        ('LYS', 'K'), ('ILE', 'I'), ('PRO', 'P'), ('THR', 'T'),
        ('PHE', 'F'), ('ASN', 'N'), ('GLY', 'G'), ('HIS', 'H'),
        ('LEU', 'L'), ('ARG', 'R'), ('TRP', 'W'), ('ALA', 'A'),
        ('VAL', 'V'), ('GLU', 'E'), ('TYR', 'Y'), ('MET', 'M'),
        # Non-canonical amino acids
        # ('MSE', 'M'), ('SOC', 'C'),
        # Canonical xNA
        ('  U', 'U'), ('  A', 'A'), ('  G', 'G'), ('  C', 'C'),
        ('  T', 'T'),
    ]

    three_to_one = dict(res_codes)
    # _records = set(['ATOM  ', 'HETATM'])
    _records = set(['ATOM  '])

    sequence = []
    read = set()
    for line in open(inpdb):
        line = line.strip()
        if line[0:6] in _records:
            resn = line[17:20]
            resi = line[22:26]
            icode = line[26]
            r_uid = (resn, resi, icode)
            if r_uid not in read:
                read.add(r_uid)
            else:
                continue
            aa_resn = three_to_one.get(resn, 'X')
            sequence.append(aa_resn)
            seq_to_res_mapping += [int(resi)]

    # return {'sequence': ''.join(sequence), 'mapping': seq_to_res_mapping}
    return ''.join(sequence)


def find_max_interface_score_from_pairs(native_pdb_pairs, pred_pdb_pairs, workdir):
    G_dockq = {}
    for native_pair_dict in native_pdb_pairs:
        native_pair = native_pair_dict['chain_ids']
        native_pdb = f"{workdir}/native_{native_pair}.pdb"
        
        chain1_pdb, chain2_pdb = native_pair_dict['pdb_str'].split(';')
        with open(native_pdb, 'w') as fw:
            for line in open(chain1_pdb):
                if line.startswith('ATOM'):
                    fw.write(line[:21] + 'A' + line[22:])
            fw.write('TER')
            for line in open(chain2_pdb):
                if line.startswith('ATOM'):
                    fw.write(line[:21] + 'B' + line[22:])
            fw.write('END')

        pair_dockq_scores = {}
        for pred_pair_dict in pred_pdb_pairs:
            pred_pair = pred_pair_dict['chain_ids']
            pred_pdb = f"{workdir}/pred_{pred_pair}.pdb"
            chain1_pdb, chain2_pdb = pred_pair_dict['pdb_str'].split(';')
            with open(pred_pdb, 'w') as fw:
                for line in open(chain1_pdb):
                    if line.startswith('ATOM'):
                        fw.write(line[:21] + 'A' + line[22:])
                fw.write('TER')
                for line in open(chain2_pdb):
                    if line.startswith('ATOM'):
                        fw.write(line[:21] + 'B' + line[22:])
                fw.write('END')
            
            pred_pair, fnat_score, dockq_score = cal_dockq_score([pred_pair, pred_pdb, native_pdb])
            pair_dockq_scores['pred_' + pred_pair] = float(dockq_score)

        G_dockq[native_pair] = pair_dockq_scores

    print(G_dockq)

    dockq_cost_matrix = []
    for native_pair in G_dockq:
        dockq_row = []
        for pair in G_dockq[native_pair]:
            dockq_row += [G_dockq[native_pair][pair]]
        dockq_cost_matrix += [dockq_row]

    dockq_cost_np = np.array(dockq_cost_matrix)
    row_ind, col_ind = linear_sum_assignment(dockq_cost_np, maximize=True)
    dockq_cost = dockq_cost_np[row_ind, col_ind].sum()

    return dockq_cost


def read_pair_info(infile):
    pdb_interact_pairs = {}
    pair = ""
    for line in open(pair_file):
        line = line.rstrip('\n')
        if len(pair) == 0:
            pair = line
        else:
            chain_pdbs = line.split(';')
            combine_seq = get_sequence(chain_pdbs[0]) + '_' + get_sequence(chain_pdbs[1])
            if combine_seq not in pdb_interact_pairs:
                pdb_interact_pairs[combine_seq] = []
            pdb_interact_pairs[combine_seq] += [dict(pdb_str=line, chain_ids=pair)]
            pair = ""
    return pdb_interact_pairs


def cal_pairwise_dockq_score_for_one_model(pdbdir, pdbname, workdir):
    pair_file = workdir + '/' + pdbname + '/pairs_info'
    pdb_interact_pairs = read_pair_info(pair_file)
    scores = []
    for other_pdb in os.listdir(pdbdir):
        if other_pdb.find('_temp') > 0:
            continue
        if other_pdb == pdbname:
            continue
        other_interact_pairs = read_pair_info(workdir + '/' + other_pdb + '/pairs_info')
        compare_dir = f"{workdir}/{pdbname}/{pdbname}_{other_pdb}"
        if not os.path.exists(compare_dir):
            os.makedirs(compare_dir)

        sum_inferface_score = 0
        other_pair_total_count = 0
        for combine_seq in other_interact_pairs:
            other_pair_total_count += len(other_interact_pairs[combine_seq])
            pair_based_interface_score = 0
            if combine_seq in pdb_interact_pairs:
                pair_based_interface_score = find_max_interface_score_from_pairs(pdb_interact_pairs[combine_seq], other_interact_pairs[combine_seq], compare_dir)
            sum_inferface_score += pair_based_interface_score

        final_interface_score = sum_inferface_score / other_pair_total_count   
        scores += [final_interface_score]
    return np.mean(np.array(scores))            


def extract_pairs(inparams):
    inpdb, outdir = inparams
    chain_pdbs = split_pdb(inpdb, outdir)
    chain_ids = list(chain_pdbs.keys())
    with open(f"{outdir}/pairs_info") as fw:
        for i in range(len(chain_ids)):
            for j in range(i+1, len(chain_ids)):
                min_dist = cal_min_distance_between_chains(chain_ids[i], chain_ids[j], chain_pdbs[chain_ids[i]]['pdb_path'], chain_pdbs[chain_ids[j]]['pdb_path'])
                if min_dist >= 6:
                    continue
                fw.write(f"{chain1}_{chain2}")
                fw.write(f"{chain_pdbs[chain_ids[i]]['pdb_path']}_{chain_pdbs[chain_ids[j]]['pdb_path']}")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    
    args = parser.parse_args()

    for target in sorted(os.listdir(args.indir)):
        
        workdir = args.outdir + '/' + target
        if not os.path.exists(workdir):
            os.makedirs(workdir)

        # Extract interaction pairs from all pdbs
        process_list = []
        for pdbfile in sorted(os.listdir(args.indir + '/' + target)):
            pdbdir = workdir + '/' + pdbfile
            if not os.path.exists(pdbdir):
                os.makedirs(pdbdir)
            process_list.append([args.indir + '/' + target + '/' + pdbfile, pdbdir])

        pool = Pool(processes=30)
        results = pool.map(extract_pairs, process_list)
        pool.close()
        pool.join()

        models = []
        scores = []
        for pdbfile in sorted(os.listdir(args.indir + '/' + target)):
            pdbdir = workdir + '/' + pdbfile
            models += [pdbfile]
            scores += [cal_pairwise_dockq_score_for_one_model(args.indir + '/' + target, pdbfile, pdbdir)]
        
        df = pd.DataFrame({'model': models, 'score': scores})
        df.to_csv(args.outdir + '/' + target + '.csv')   

        


           



            
            





    

    
