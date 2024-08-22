import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr

def read_chain_contents(inpdb):
    chain_contents = {}
    for line in open(inpdb, 'r').readlines():
        if not line.startswith('ATOM'):
            continue
        chain_name = line[21]
        if chain_name not in chain_contents:
            chain_contents[chain_name] = []
        chain_contents[chain_name] += [line]
    return chain_contents

def reindex_pdb(inpdb, outpdb):
    resCounter = 0
    atomCounter = 0
    prevrNum = "XX"
    prevchain = "XX"
    contents = []
    for line in open(inpdb):
        if not line.startswith('ATOM'):
            continue
        rnum = line[22:27]
        chain = line[21:22]
        add_ter = False
        if prevchain != chain:
            if prevchain != "XX":
                contents += ["TER\n"]
            prevchain = chain
            resCounter = 1
            atomCounter = 0
            prevrNum = rnum
        elif prevrNum != rnum:
            prevrNum = rnum
            resCounter += 1
        atomCounter += 1
        rnum_string = "{:>4}".format(resCounter)
        anum_string = "{:>5}".format(atomCounter)
        if chain == "0":
            row_chain = "a"
        else:
            row_chain = chain
        row = f"{line[:6]}{anum_string}{line[11:21]}{row_chain}{rnum_string}{line[26:]}"
        contents += [row]
    contents += ["TER\n"]
    with open(outpdb, 'w') as fw:
        fw.writelines("".join(contents))
        fw.write('END')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--inpdbs', type=str, required=True)
    parser.add_argument('--outpdb', type=str, required=True)
    args = parser.parse_args()

    with open(args.outpdb + '.temp', 'w') as fw:
        pdb_contents = []
        for inpdb in args.inpdbs.split(','):
            chain_contents = read_chain_contents(inpdb)
            pdb_contents += [chain_contents]
        
        chain_ids = set()
        for pdb_content in pdb_contents:
            for chain_id in pdb_content:
                chain_ids.add(chain_id)

        chain_ids = sorted(list(chain_ids))
        print(chain_ids)
        for chain_id in sorted(chain_ids):
            contents = []
            for pdb_content in pdb_contents:
                if chain_id in pdb_content:
                    contents += pdb_content[chain_id]
            fw.writelines(''.join(contents))
            fw.write('TER\n')
        fw.write('END')

    reindex_pdb(args.outpdb + '.temp', args.outpdb)
