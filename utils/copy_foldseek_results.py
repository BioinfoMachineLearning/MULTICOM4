import os, sys, argparse, time, pathlib
from multiprocessing import Pool
from multicom4.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from multicom4.monomer_alignment_generation.alignment import write_fasta
from multicom4.common.protein import read_qa_txt_as_df, parse_fasta, complete_result, make_chain_id_map
from multicom4.multimer_structure_refinement import iterative_refine_pipeline_multimer
from multicom4.monomer_structure_refinement import iterative_refine_pipeline
from multicom4.common.pipeline import run_monomer_msa_pipeline, run_monomer_template_search_pipeline, \
    run_monomer_structure_generation_pipeline_v2, run_monomer_evaluation_pipeline, \
    run_monomer_msas_concatenation_pipeline, run_monomer_templates_concatenation_pipeline, \
    run_multimer_structure_generation_pipeline_v2, \
    run_multimer_structure_generation_pipeline_foldseek, \
    run_multimer_evaluation_pipeline, run_monomer_msa_pipeline_img, foldseek_iterative_monomer_input, \
    copy_same_sequence_msas
from absl import flags
from absl import app
import copy
import pandas as pd
from multicom4.common.config import *

flags.DEFINE_string('indir', None, 'option file')
flags.DEFINE_string('stoi', None, 'Path to monomer fasta')
FLAGS = flags.FLAGS


def main(argv):

    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    src_methods = ['folds_iter_not', 'folds_iter_esm_not']

    heteromer_run_methods = ['folds_iter', 'folds_iter_nop', 'folds_iter_notp', 'folds_iter_esm', 'folds_iter_esm_nop', 'folds_iter_esm_notp']

    homomer_run_methods = ['folds_iter', 'folds_iter_o', 'folds_iter_o_not',
                            'folds_iter_esm', 'folds_iter_esm_o', 'folds_iter_esm_o_not']

    run_methods = homomer_run_methods if FLAGS.stoi == "homomer" else heteromer_run_methods

    for run_method in run_methods:
        src_method = 'folds_iter_not'
        if run_method.find('esm') >= 0:
            src_method = 'folds_iter_esm_not'

        for i in range(1, 3):
            srcdir = os.path.join(FLAGS.indir, f'{src_method}_{i}', 'prepare')

            trgdir = os.path.join(FLAGS.indir, f'{run_method}_{i}')
            os.makedirs(trgdir, exist_ok=True)

            trgdir = os.path.join(trgdir, 'prepare')
            if os.path.exists(trgdir):
                os.system(f"rm -rf {trgdir}")

            os.system(f"cp -r {srcdir} {trgdir}")
    

if __name__ == '__main__':
    flags.mark_flags_as_required([
        'indir',
        'stoi',
    ])
    app.run(main)
