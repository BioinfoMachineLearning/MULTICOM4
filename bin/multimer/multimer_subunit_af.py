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

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to multimer fasta')
flags.DEFINE_string('output_dir', None, 'Output directory')
flags.DEFINE_boolean('run_img', False, 'Whether to use IMG alignment to generate models')
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

    processed_seuqences = {}
    monomer_run_methods = ['default', 'deepmsa_DeepJGI', 'deepmsa_DeepJGI_hms', 'deepmsa_aMSA', 'deepmsa_dMSA', 'deepmsa_dMSA_hhb',
                            'deepmsa_dMSA_hms', 'deepmsa_dMSA_jac', 'deepmsa_q3JGI', 'deepmsa_q3JGI_hms', 'deepmsa_q4JGI',
                            'deepmsa_q4JGI_hms', 'deepmsa_qMSA', 'deepmsa_qMSA_hh3', 'deepmsa_qMSA_hhb', 'deepmsa_qMSA_hms',
                            'deepmsa_qMSA_jac', 'def_drop_nos', 'def_drop_s', 
                            'def_esm_msa', 'def_esm_msa_ckpt5', 
                            'def_notemp',
                            'def_notemp_drop_nos', 'def_notemp_drop_s', 'def_ptm_not_drop_s', 'def_ptm_drop_nos', 'def_ptm_drop_s',
                            'def_ptm_not_drop_nos', 'def_ptm_notemp', 'default_pdb70_new', 'default_ptm',
                            'default_seq_temp', 'dhr', 'ori_seq_temp', 'original']

    for chain_id in chain_id_map:
        monomer_id = chain_id
        monomer_sequence = chain_id_map[chain_id].sequence
        monomer_fasta = os.path.join(FLAGS.output_dir, monomer_id + ".fasta")

        if monomer_sequence not in processed_seuqences:
            with open(monomer_fasta, "w") as fw:
                write_fasta({chain_id: monomer_sequence}, fw)
            N1_monomer_outdir = os.path.join(N1_outdir, monomer_id)
            N1_monomer_outdir_img = os.path.join(N1_outdir_img, monomer_id)

            N2_monomer_outdir = os.path.join(N2_outdir, monomer_id)

            N3_monomer_outdir = os.path.join(N3_outdir, monomer_id)
            makedir_if_not_exists(N3_monomer_outdir)
            if not run_monomer_structure_generation_pipeline_v2(params=params,
                                                                targetname=targetname,
                                                                fasta_path=monomer_fasta,
                                                                alndir=N1_monomer_outdir,
                                                                img_alndir=N1_monomer_outdir_img,
                                                                templatedir=N2_monomer_outdir,
                                                                outdir=N3_monomer_outdir,
                                                                run_methods=monomer_run_methods,
                                                                run_script=True,
                                                                is_subunit=True):
                print(f"Program failed in step 3: monomer {monomer_id} structure generation")

            processed_seuqences[monomer_sequence] = monomer_id


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
