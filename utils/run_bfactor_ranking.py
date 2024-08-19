import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
from biopandas.pdb import PandasPdb

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--inpdb', type=str, required=True)
    args = parser.parse_args()

    ppdb = PandasPdb()
    ppdb.read_pdb(args.inpdb)
    plddt = ppdb.df['ATOM'][ppdb.df['ATOM']['atom_name'] == 'CA']['b_factor']
    plddt = plddt.to_numpy().astype(np.float32)
    print(np.mean(np.array(plddt)))

