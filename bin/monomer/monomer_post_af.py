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
from multicom4.common.config import *

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to monomer fasta')
flags.DEFINE_string('output_dir', None, 'Output directory')
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
    makedir_if_not_exists(N1_outdir)

    N1_outdir_img = os.path.join(outdir, 'N1_monomer_alignments_generation_img')

    print("#################################################################################################")
    print(f"1. Start to generate alignments for monomers")

    #result = run_monomer_msa_pipeline(fasta=FLAGS.fasta_path, outdir=N1_outdir, params=params, only_monomer=True)

    #if result is None:
    #    raise RuntimeError('The monomer alignment generation has failed!')

    print("#################################################################################################")
    print("2. Start to generate template for monomer")

    N2_outdir = os.path.join(outdir, 'N2_monomer_template_search')

    makedir_if_not_exists(N2_outdir)

    template_file = run_monomer_template_search_pipeline(params=params, targetname=targetname, sequence=sequence,
                                                         a3m=os.path.join(N1_outdir, targetname+"_uniref90.sto"), outdir=N2_outdir)

    if template_file is None:
        raise RuntimeError("Program failed in step 2: monomer template search")

    print("The generation for monomer template has finished!")

    print("3. Start to generate tertiary structure for monomers using alphafold")
    N3_outdir = os.path.join(outdir, 'N3_monomer_structure_generation')
    makedir_if_not_exists(N3_outdir)
    run_methods = [ 'def_esm_msa', 'def_esm_msa_ckpt5']

    if not run_monomer_structure_generation_pipeline_v2(params=params,
                                                        targetname=targetname,
                                                        fasta_path=FLAGS.fasta_path,
                                                        alndir=N1_outdir, 
                                                        img_alndir=N1_outdir_img,
                                                        templatedir=N2_outdir, 
                                                        outdir=N3_outdir,
                                                        run_methods=run_methods,
                                                        run_script=True):
        print("Program failed in step 3: monomer structure generation")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
