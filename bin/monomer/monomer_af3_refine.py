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
import json

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to monomer fasta')
flags.DEFINE_string('output_dir', None, 'Output directory')
flags.DEFINE_string('config_name', 'af3_refine', 'single predictor name')
flags.DEFINE_string('inpdb', None, 'Whether to use IMG alignment to generate models')
FLAGS = flags.FLAGS


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
    N1_outdir_img = os.path.join(outdir, 'N1_monomer_alignments_generation_img')
    N2_outdir = os.path.join(outdir, 'N2_monomer_template_search')

    print("#################################################################################################")

    print("#################################################################################################")

    print("3. Start to generate tertiary structure for monomers using alphafold")
    N3_outdir = os.path.join(outdir, 'N3_monomer_structure_generation')
    makedir_if_not_exists(N3_outdir)
    
    refinement_inputs= []
    refine_input = iterative_refine_pipeline.refinement_input(fasta_path=FLAGS.fasta_path,
                                                              pdb_path=FLAGS.inpdb,
                                                              pkl_path='',
                                                              msa_path='')
    refinement_inputs += [refine_input]

    refine_dir = os.path.join(N3_outdir, FLAGS.config_name, 'workdir')
    makedir_if_not_exists(refine_dir)
    pipeline = iterative_refine_pipeline.Monomer_iterative_refinement_pipeline_server(params=params, config_name=FLAGS.config_name)
    pipeline.search(refinement_inputs=refinement_inputs, outdir=refine_dir,
                    uniref90_sto=os.path.join(N1_outdir, targetname + '_uniref90.sto'))
    

if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir',
    ])
    app.run(main)