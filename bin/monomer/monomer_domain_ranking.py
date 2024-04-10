import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from multicom4.monomer_templates_concatenation import parsers
from multicom4.monomer_alignment_generation.alignment import *
from multicom4.monomer_structure_evaluation.bfactor_ranking import Bfactor_qa

import os, sys, argparse, time, copy, pathlib
from multiprocessing import Pool
from multicom4.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from multicom4.monomer_structure_refinement import iterative_refine_pipeline
from multicom4.common.protein import read_qa_txt_as_df, complete_result
from multicom4.common.pipeline import run_monomer_msa_pipeline, run_monomer_template_search_pipeline, \
    run_monomer_structure_generation_pipeline_v2, run_monomer_evaluation_pipeline, \
    run_monomer_msa_pipeline_img
import pandas as pd
from absl import flags
from absl import app

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to monomer fasta')
flags.DEFINE_string('output_dir', None, 'Output directory')
flags.DEFINE_string('domain_info', None, 'Output directory')
FLAGS = flags.FLAGS

def find_idx_in_domain_ranges(residue_idx, domain_starts, domain_ends):
    for start, end in zip(domain_starts, domain_ends):
        if start <= residue_idx <= end:
            return True
    return False

def main(argv):

    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    os.environ['TF_FORCE_UNIFIED_MEMORY'] = '1'
    os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '4.0'

    FLAGS.fasta_path = os.path.abspath(FLAGS.fasta_path)
    FLAGS.output_dir = os.path.abspath(FLAGS.output_dir)

    check_file(FLAGS.option_file)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    check_file(FLAGS.fasta_path)

    targetname = pathlib.Path(FLAGS.fasta_path).stem
    sequence = open(FLAGS.fasta_path).readlines()[1].rstrip('\n').strip()

    outdir = FLAGS.output_dir

    makedir_if_not_exists(outdir)

    N1_outdir = os.path.join(outdir, 'N1_monomer_alignments_generation')
    makedir_if_not_exists(N1_outdir)

    contact_map_file = os.path.join(N1_outdir, 'dncon4', f'{targetname}.dncon2.rr')
    if not os.path.exists(contact_map_file):
        raise Exception("The contact map file hasn't been generated!")

    dist_map_file = os.path.join(N1_outdir, 'deepdist', f'{targetname}.txt')
    if not os.path.exists(dist_map_file):
        raise Exception("The distance map file hasn't been generated!")

    print("#################################################################################################")
    print(f"1. Start to generate alignments for monomers")
    
    print("#################################################################################################")
    print("2. Start to generate template for monomer")

    N2_outdir = os.path.join(outdir, 'N2_monomer_template_search')

    print("The generation for monomer template has finished!")

    print("#################################################################################################")

    print("#################################################################################################")

    print("3. Start to generate tertiary structure for monomers using alphafold")

    N3_outdir = os.path.join(outdir, 'N3_monomer_structure_generation')

    print("#################################################################################################")

    print("#################################################################################################")

    print("4. Start to evaluate monomer models")

    N4_outdir = os.path.join(outdir, 'N4_monomer_domain_structure_evaluation')

    makedir_if_not_exists(N4_outdir)
    
    all_domain_ranges = []

    for line in open(FLAGS.domain_info):
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
    
    full_pdbdir = os.path.join(N4_outdir, 'full')
    os.makedirs(full_pdbdir, exist_ok=True)

    to_be_ranked_dirs = [full_pdbdir]

    bfactorqa = Bfactor_qa()
    model_count = 5
    for method in os.listdir(N3_outdir):

        ranking_json_file = os.path.join(N3_outdir, method, "ranking_debug.json")
        if not os.path.exists(ranking_json_file):
            continue

        ranking_json = json.loads(open(ranking_json_file).read())
        for i in range(model_count):
            pdbname = f"{method}_{i}"
            ranked_pdb = os.path.join(N3_outdir, method, f"ranked_{i}.pdb")
            trg_pdb = os.path.join(N3_outdir, f"{pdbname}.pdb")
            os.system(f"cp {ranked_pdb} {trg_pdb}")
                
            contents = open(ranked_pdb).readlines()
            for domain_idx, (domain_starts, domain_ends) in enumerate(all_domain_ranges):
                domain_pdbdir = os.path.join(N4_outdir, 'domain' + str(domain_idx))
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
                domain_ranked_pdb = os.path.join(domain_pdbdir, f"{pdbname}.pdb")
                open(domain_ranked_pdb, 'w').writelines(''.join(domain_contents))
    
    for to_be_ranked_dir in to_be_ranked_dirs:
        ranking_file = to_be_ranked_dir + '_plddt.csv'
        bfactor_ranking = bfactorqa.run(input_dir=to_be_ranked_dir)
        bfactor_ranking.to_csv(ranking_file)
    
    # may need to add average domain ranking later

    if result is None:
        raise RuntimeError("Program failed in step 4: monomer model evaluation")

    print("The evaluation for monomer models has been finished!")

    print("#################################################################################################")

    print("#################################################################################################")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir',
    ])
    app.run(main)

    
    