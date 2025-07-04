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
    run_multimer_evaluation_pipeline, foldseek_iterative_monomer_input, \
    copy_same_sequence_msas

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
    targetname = os.path.basename(FLAGS.fasta_path).replace('.fasta', '')
    
    check_file(FLAGS.option_file)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    check_file(FLAGS.fasta_path)

    makedir_if_not_exists(FLAGS.output_dir)

    N1_outdir = os.path.join(FLAGS.output_dir, 'N1_monomer_alignments_generation')
    N2_outdir = os.path.join(FLAGS.output_dir, 'N2_monomer_template_search')
    N3_outdir = os.path.join(FLAGS.output_dir, 'N3_monomer_structure_generation')
    
    print("#################################################################################################")

    print("#################################################################################################")
    print("1-3. Start to generate monomer models")

    makedir_if_not_exists(N1_outdir)

    with open(FLAGS.fasta_path) as f:
        input_fasta_str = f.read()
    input_seqs, input_descs = parse_fasta(fasta_string=input_fasta_str)
    chain_id_map, chain_id_seq_map = make_chain_id_map(sequences=input_seqs,
                                                       descriptions=input_descs)

    processed_seuqences = {}

    for chain_id in chain_id_map:
        monomer_id = chain_id
        monomer_sequence = chain_id_map[chain_id].sequence
        monomer_fasta = os.path.join(FLAGS.output_dir, monomer_id + ".fasta")

        if monomer_sequence not in processed_seuqences:
            with open(monomer_fasta, "w") as fw:
                write_fasta({chain_id: monomer_sequence}, fw)
            N1_monomer_outdir = os.path.join(N1_outdir, monomer_id)
            makedir_if_not_exists(N1_monomer_outdir)
            result = run_monomer_msa_pipeline(monomer_fasta, N1_monomer_outdir, params)
            if result is None:
                raise RuntimeError(f"Program failed in step 1: monomer {monomer_id} alignment generation")

            N2_monomer_outdir = os.path.join(N2_outdir, monomer_id)
            makedir_if_not_exists(N2_monomer_outdir)
            template_file = run_monomer_template_search_pipeline(targetname=monomer_id, sequence=monomer_id,
                                                                 a3m=os.path.join(N1_monomer_outdir, monomer_id + "_uniref90.sto"),
                                                                 outdir=N2_monomer_outdir, params=params)
            if template_file is None:
                raise RuntimeError(f"Program failed in step 2: monomer {monomer_id} template search")

            processed_seuqences[monomer_sequence] = monomer_id



if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
