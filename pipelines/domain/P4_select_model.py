import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from multicom4.monomer_templates_concatenation import parsers
from multicom4.monomer_alignment_generation.alignment import *
from multicom4.monomer_structure_evaluation.bfactor_ranking import Bfactor_qa

def find_idx_in_domain_ranges(residue_idx, domain_starts, domain_ends):
    for start, end in zip(domain_starts, domain_ends):
        if start <= residue_idx <= end:
            return True
    return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #parser.add_argument('--fasta_path', type=str, required=True)
    parser.add_argument('--domain_info', type=str, required=True)
    parser.add_argument('--domain_run_dir', type=str, required=True)
    parser.add_argument('--outputdir', type=str, required=True)
    args = parser.parse_args()

    all_domain_ranges = []

    for line in open(args.domain_info):
        line = line.rstrip('\n')
        # domain 1: 140-317 Normal
        # domain 1: 140-317
        domain_name, domain_range = line.split(':')
        domain_name = domain_name.replace(' ', '')
        domain_range = domain_range.split()[0]

        domain_sequence = ""
        starts, ends = [], []

        for domain_range_str in domain_range.split(','):
            start, end = domain_range_str.split('-')
            # list index start from 0
            starts += [int(start)]
            ends += [int(end)]

        all_domain_ranges += [(starts, ends)]

    os.makedirs(args.outputdir, exist_ok=True)

    method_dirs = os.listdir(args.domain_run_dir)

    pdbdir = os.path.join(args.outputdir, 'pdb')
    os.makedirs(pdbdir, exist_ok=True)
    
    full_pdbdir = os.path.join(pdbdir, 'full')
    os.makedirs(full_pdbdir, exist_ok=True)

    to_be_ranked_dirs = [full_pdbdir]

    bfactorqa = Bfactor_qa()
    for method_dir in method_dirs:
        resultdir = os.path.join(args.domain_run_dir, method_dir, 'alphafold')
        for ranked_pdb in os.listdir(resultdir):
            if ranked_pdb.find('ranked_') < 0:
                continue
            
            full_ranked_pdb = os.path.join(full_pdbdir, method_dir + '_' + ranked_pdb)
            os.system(f"cp {resultdir}/{ranked_pdb} {full_ranked_pdb}")

            contents = open(full_ranked_pdb).readlines()

            for domain_idx, (domain_starts, domain_ends) in enumerate(all_domain_ranges):
                domain_pdbdir = os.path.join(pdbdir, 'domain' + str(domain_idx))
                os.makedirs(domain_pdbdir, exist_ok=True)
                if domain_pdbdir not in to_be_ranked_dirs:
                    to_be_ranked_dirs += [domain_pdbdir]

                # may need to reorder the residue index later
                domain_contents = []
                for line in contents:
                    if not line.startswith('ATOM'):
                        continue
                    residue_num = int(line[22:26].strip())
                    if find_idx_in_domain_ranges(residue_num, domain_starts, domain_ends):
                        domain_contents += [line]
                domain_ranked_pdb = os.path.join(domain_pdbdir, method_dir + '_' + ranked_pdb)
                open(domain_ranked_pdb, 'w').writelines(''.join(domain_contents))
    
    for to_be_ranked_dir in to_be_ranked_dirs:
        ranking_file = to_be_ranked_dir + '_plddt.csv'
        bfactor_ranking = bfactorqa.run(input_dir=to_be_ranked_dir)
        bfactor_ranking.to_csv(ranking_file)
    
    # may need to add average domain ranking later
    