import os, sys, argparse, time
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
    copy_same_sequence_msas, select_final_monomer_models_from_complex

from absl import flags
from absl import app
import copy
import pandas as pd

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to multimer fasta')
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

    makedir_if_not_exists(FLAGS.output_dir)

    N1_outdir = os.path.join(FLAGS.output_dir, 'N1_monomer_alignments_generation')
    N1_outdir_img = os.path.join(FLAGS.output_dir, 'N1_monomer_alignments_generation_img') 
    N2_outdir = os.path.join(FLAGS.output_dir, 'N2_monomer_template_search')
    N3_outdir = os.path.join(FLAGS.output_dir, 'N3_monomer_structure_generation')
    img_msas = {}

    print("#################################################################################################")

    print("#################################################################################################")
    print("1-3. Start to generate monomer models")

    makedir_if_not_exists(N1_outdir)

    with open(FLAGS.fasta_path) as f:
        input_fasta_str = f.read()
    input_seqs, input_descs = parse_fasta(fasta_string=input_fasta_str)
    chain_id_map, chain_id_seq_map = make_chain_id_map(sequences=input_seqs,
                                                       descriptions=input_descs)

    print("#################################################################################################")

    print("#################################################################################################")
    print("4. Start to generate complex alignments")

    N4_outdir = os.path.join(FLAGS.output_dir, 'N4_monomer_alignments_concatenation')
    
    print("#################################################################################################")

    print("#################################################################################################")

    print("5. Start to search complex templates based on monomer structures")

    N5_outdir = os.path.join(FLAGS.output_dir, 'N5_monomer_templates_concatenation')

    print("#################################################################################################")

    print("#################################################################################################")

    print("6. Start to generate complex multimer structures")

    N6_outdir = os.path.join(FLAGS.output_dir, 'N6_multimer_structure_generation')

    run_methods = ['def_mul_refine', 'af3_refine']
    
    stoichiometry = "homomers" if len(set(input_seqs)) == 1 else 'heteromers'
    for run_method in run_methods:
        if os.path.exists(os.path.join(N6_outdir, run_method)):
            if not complete_result(os.path.join(N6_outdir, run_method), 5):
                refine_dir = os.path.join(N6_outdir, run_method, 'workdir')
                final_dir = os.path.join(N6_outdir, run_method, 'finaldir')
                os.makedirs(final_dir, exist_ok=True)
                pipeline = iterative_refine_pipeline_multimer.Multimer_refinement_model_selection(params=params, config_name=run_method, stoichiometry=stoichiometry)
                pipeline.select_v1(indir=refine_dir, outdir=final_dir)
                pipeline.make_predictor_results(final_dir, os.path.join(N6_outdir, run_method))

    print("Multimer structure generation has been finished!")

    print("#################################################################################################")

    print("#################################################################################################")

    print("9. Start to evaluate multimer models")
    run_methods=["alphafold", "bfactor", 'multieva', 'gate', 'gate_noaf']
    N7_outdir = os.path.join(FLAGS.output_dir, 'N7_multimer_structure_evaluation')
    multimer_qa_result = run_multimer_evaluation_pipeline(fasta_path=FLAGS.fasta_path,
                                                          params=params, run_methods=run_methods,
                                                          chain_id_map=chain_id_map,
                                                          indir=N6_outdir, outdir=N7_outdir,
                                                          is_human=True)

    print("#################################################################################################")

    # N8_outdir = os.path.join(FLAGS.output_dir, 'N8_monomer_structure_evaluation')
    
    # select_final_monomer_models_from_complex(chain_id_map=chain_id_map, 
    #                                          multimer_qa_result_dir=N7_outdir, 
    #                                          outputdir=N8_outdir)

    print("#################################################################################################")

if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
