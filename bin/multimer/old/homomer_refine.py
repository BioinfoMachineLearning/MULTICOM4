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
    copy_same_sequence_msas

from absl import flags
from absl import app
import copy
import pandas as pd
import json 

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to multimer fasta')
flags.DEFINE_string('output_dir', None, 'Output directory')
flags.DEFINE_string('config_name', None, 'Whether to use IMG alignment to generate models')
flags.DEFINE_integer('idx', None, 'Whether to use IMG alignment to generate models')

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

    method_out_dir = os.path.join(N6_outdir, FLAGS.config_name) 
    # refine default top-ranked models
    refinement_inputs = []
    default_workdir = os.path.join(N6_outdir, 'default_multimer')
    ranking_json_file = os.path.join(default_workdir, "ranking_debug.json")
    if not os.path.exists(ranking_json_file):
        raise Exception(f"Haven't generated default_multimer models!")

    ranking_json = json.loads(open(ranking_json_file).read())

    i = FLAGS.idx
    pdb_path = os.path.join(default_workdir, f"ranked_{i}.pdb")
    model_name = list(ranking_json["order"])[i]
    pkl_path = os.path.join(default_workdir, f"result_{model_name}.pkl")
    msa_paths = {}
    for chain_id in chain_id_map:
        monomer_msa = os.path.join(default_workdir, 'msas', chain_id, "monomer_final.a3m")
        paired_msa = os.path.join(default_workdir, 'msas', f"{chain_id}.paired.a3m")
        msa_paths[chain_id] = dict(paired_msa=paired_msa,
                                    monomer_msa=monomer_msa)

    refine_input = iterative_refine_pipeline_multimer.refinement_input_multimer(chain_id_map=chain_id_map,
                                                                                fasta_path=FLAGS.fasta_path,
                                                                                pdb_path=pdb_path,
                                                                                pkl_path=pkl_path,
                                                                                msa_paths=msa_paths)
    refinement_inputs += [refine_input]

    refine_dir = os.path.join(method_out_dir, 'workdir')
    makedir_if_not_exists(refine_dir)

    pipeline = iterative_refine_pipeline_multimer.Multimer_iterative_refinement_pipeline_server(params=params, config_name=FLAGS.config_name)
    pipeline.search(refinement_inputs=refinement_inputs, outdir=refine_dir, stoichiometry="homomer")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
