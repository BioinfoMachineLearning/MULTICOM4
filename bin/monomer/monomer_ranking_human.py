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

    contact_map_file = os.path.join(N1_outdir, 'dncon4', f'{targetname}.dncon2.rr')
    #if not os.path.exists(contact_map_file):
    #    raise Exception(f"The contact map file hasn't been generated: {contact_map_file}")

    dist_map_file = os.path.join(N1_outdir, 'deepdist', f'{targetname}.txt')
    #if not os.path.exists(dist_map_file):
    #    raise Exception(f"The distance map file hasn't been generated: {dist_map_file}")

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

    # run_methods = ['foldseek_refine', 'foldseek_refine_esm', 'foldseek_refine_esm_h']
    run_methods = []
    for run_method in run_methods:
        if not complete_result(os.path.join(N3_outdir, run_method), 5):
            refine_dir = os.path.join(N3_outdir, run_method, 'workdir')
            final_dir = os.path.join(N3_outdir, run_method, 'finaldir')
            pipeline = iterative_refine_pipeline.Monomer_refinement_model_selection(params=params, config_name=run_method)
            pipeline.select_v1(indir=refine_dir, outdir=final_dir)
            pipeline.make_predictor_results(final_dir, os.path.join(N3_outdir, run_method))

    print("The prediction for monomers has finished!")

    print("#################################################################################################")

    print("#################################################################################################")

    print("4. Start to evaluate monomer models")

    N4_outdir = os.path.join(outdir, 'N4_monomer_structure_evaluation')

    makedir_if_not_exists(N4_outdir)

    run_methods = None
    if os.path.exists(params['slurm_script_template']):
        run_methods = ["alphafold", "bfactor", "apollo"]

    result = run_monomer_evaluation_pipeline(params=params, targetname=targetname, fasta_file=FLAGS.fasta_path,
                                             input_monomer_dir=N3_outdir, outputdir=N4_outdir,
                                             generate_final_models=True, contact_map_file=contact_map_file,
                                             dist_map_file=dist_map_file, run_methods=run_methods, is_human=True)

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
