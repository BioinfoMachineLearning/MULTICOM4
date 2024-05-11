import os, sys, argparse, time, copy, pathlib
from multiprocessing import Pool
from multicom4.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from multicom4.common.protein import read_qa_txt_as_df, complete_result
from multicom4.common.pipeline import run_monomer_msa_pipeline, run_monomer_template_search_pipeline
import pandas as pd
from absl import flags
from absl import app

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to multimer fastas')
flags.DEFINE_string('output_dir', None, 'Output directory')
FLAGS = flags.FLAGS

def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    check_file(FLAGS.option_file)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    params['colabfold_search_program'] = ""
    params['colabfold_split_msas_program'] = ""
    params['colabfold_databases'] = ""
    params['mmseq_program'] = ""

    params['deepmsa2_path'] = ""
    params['JGIclust_database'] = ""
    params['metaclust_database'] = ""
    
    params['DHR_program_path'] = ""
    params['DHR_database_path'] = ""

    # targetname = pathlib.Path(FLAGS.fasta_path).stem
    targetname = None
    sequence = None
    for line in open(FLAGS.fasta_path):
        line = line.rstrip('\n').strip()
        if line.startswith('>'):
            targetname = line[1:].split()[0]
            # if targetname_in_fasta != targetname:
            #     print("Warning: fasta file name doesn't match with fasta content!")
        else:
            sequence = line

    makedir_if_not_exists(FLAGS.output_dir)

    N1_outdir = FLAGS.output_dir + '/N1_monomer_alignments_generation'
    makedir_if_not_exists(N1_outdir)

    print("#################################################################################################")
    print(f"1. Start to generate alignments for monomers")

    params['uniclust_db'] = []
    params['colabfold_databases'] = []
    params['JGIclust_database'] = []
    params['DHR_database_path'] = []

    result = run_monomer_msa_pipeline(fasta=FLAGS.fasta_path, outdir=N1_outdir, 
                                      params=params, only_monomer=True, run_auxi_output=False)

    if result is None:
        raise RuntimeError('The monomer alignment generation has failed!')


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
