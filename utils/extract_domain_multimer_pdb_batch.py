import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr

def split_pdb(complex_pdb, outdir):
    pre_chain = None
    i = 0
    for line in open(complex_pdb, 'r').readlines():
        if not line.startswith('ATOM'):
            continue
        chain_name = line[21]
        if pre_chain is None:
            pre_chain = chain_name
            fw = open(outdir + '/' + chain_name + '.pdb', 'w')
            fw.write(line[:21] + ' ' + line[22:])
        elif chain_name == pre_chain:
            fw.write(line[:21] + ' ' + line[22:])
        else:
            fw.close()
            i = i + 1
            fw = open(outdir + '/' + chain_name + '.pdb', 'w')
            fw.write(line[:21] + ' ' + line[22:])
            pre_chain = chain_name
    fw.close()


def combine_pdb(indir, outpdb):
    with open(outpdb, 'w') as fw:
        for chain_pdb in sorted(os.listdir(indir)):
            if chain_pdb.find('reindex') > 0:
                for line in open(indir + '/' + chain_pdb):
                    chain_id = chain_pdb[0]
                    if line.startswith('ATOM'):
                        fw.write(line[:21] + chain_id + line[22:])
                fw.write("TER\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str, required=True)
    parser.add_argument('--workdir', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    parser.add_argument('--domain_infos', type=str, required=True)
    parser.add_argument('--start_idx', type=str, required=True)
    args = parser.parse_args()

    extract_script = '/bmlfast/bml_casp16/MULTICOM4/utils/extract_domain.pl'
    reindex_script = '/bmlfast/bml_casp16/MULTICOM4/utils/domain_extract/reindex_pdb_index.pl'

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.workdir, exist_ok=True)

    for inpdb in os.listdir(args.indir):
        workdir = os.path.join(args.workdir, inpdb)
        os.makedirs(workdir, exist_ok=True)
        split_pdb(os.path.join(args.indir, inpdb), workdir)
        for domain_info, start_idx in zip(args.domain_infos.split(','), args.start_idx.split(',')):
            start_idx = int(start_idx)
            chain, chain_start, chain_end = domain_info.split('_')
            os.system(f"perl {extract_script} {workdir}/{chain}.pdb {workdir}/{chain}.domain.pdb "
                        f"{chain_start} {chain_end}")
            os.system(f"perl {reindex_script} {workdir}/{chain}.domain.pdb {workdir}/{chain}.reindex.pdb {start_idx-1}")
        combine_pdb(workdir, workdir + '/final.pdb')
        os.system(f"cp {workdir}/final.pdb {args.outdir}/{inpdb}")
