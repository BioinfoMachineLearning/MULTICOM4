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


def get_interact_pos_between_chains(chain1, chain2, chain1_pdb, chain2_pdb):
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


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--inpdb', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)

    args = parser.parse_args()

    chain_pdbs = split_pdb(args.inpdb, args.outdir)



           



            
            





    

    
