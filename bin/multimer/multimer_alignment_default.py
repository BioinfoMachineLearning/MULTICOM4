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
flags.DEFINE_boolean('run_img', True, 'Whether to use IMG alignment to generate models')
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
    monomer_run_methods = ['default', 'default_seq_temp','def_drop_s','def_drop_nos',
                            'def_notemp', 'def_notemp_drop_s', 'def_notemp_drop_nos',
                            'original', 'ori_seq_temp', #'colabfold', 'colab_seq_temp',
                            #'img', 'img_seq_temp', 
                            'dhr', 
                            'deepmsa_dMSA_hhb', 'deepmsa_dMSA_jac', 'deepmsa_dMSA_hms',
                            'deepmsa_dMSA', 'deepmsa_qMSA', 'deepmsa_aMSA', 'deepmsa_qMSA_hhb',
                            'deepmsa_qMSA_jac', 'deepmsa_qMSA_hh3', 'deepmsa_qMSA_hms',
                            'deepmsa_DeepJGI_hms', 'deepmsa_DeepJGI', 'deepmsa_q3JGI', 
                            'deepmsa_q4JGI', 'deepmsa_q3JGI_hms', 'deepmsa_q4JGI_hms',
                            'paddle-helix', 'esmfold', 'deepfold', 'megafold', 'def_esm_msa']

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

            if FLAGS.run_img:
                N1_monomer_outdir_img = os.path.join(N1_outdir_img, monomer_id)
                makedir_if_not_exists(N1_monomer_outdir_img)
                img_msas[chain_id] = run_monomer_msa_pipeline_img(params=params,
                                                                fasta=monomer_fasta,
                                                                outdir=N1_monomer_outdir_img)

            N2_monomer_outdir = os.path.join(N2_outdir, monomer_id)
            makedir_if_not_exists(N2_monomer_outdir)
            template_file = run_monomer_template_search_pipeline(targetname=monomer_id, sequence=monomer_id,
                                                                 a3m=os.path.join(N1_monomer_outdir, monomer_id + "_uniref90.sto"),
                                                                 outdir=N2_monomer_outdir, params=params)
            if template_file is None:
                raise RuntimeError(f"Program failed in step 2: monomer {monomer_id} template search")

            N3_monomer_outdir = os.path.join(N3_outdir, monomer_id)
            makedir_if_not_exists(N3_monomer_outdir)
            if not run_monomer_structure_generation_pipeline_v2(params=params,
                                                                targetname=targetname,
                                                                fasta_path=monomer_fasta,
                                                                alndir=N1_monomer_outdir,
                                                                img_alndir=N1_monomer_outdir_img,
                                                                templatedir=N2_monomer_outdir,
                                                                outdir=N3_monomer_outdir,
                                                                run_methods=monomer_run_methods):
                print(f"Program failed in step 3: monomer {monomer_id} structure generation")

            processed_seuqences[monomer_sequence] = monomer_id
        else:
            with open(monomer_fasta, "w") as fw:
                write_fasta({chain_id: monomer_sequence}, fw)
            N1_monomer_outdir = os.path.join(N1_outdir, monomer_id)
            makedir_if_not_exists(N1_monomer_outdir)

            copy_same_sequence_msas(srcdir=os.path.join(N1_outdir, processed_seuqences[monomer_sequence]),
                                    trgdir=N1_monomer_outdir,
                                    srcname=processed_seuqences[monomer_sequence],
                                    trgname=monomer_id)

            N1_monomer_deepmsa_outdir = os.path.join(N1_outdir, monomer_id, 'DeepMSA2_a3m', 'finalMSAs')
            makedir_if_not_exists(N1_monomer_deepmsa_outdir)

            copy_same_sequence_msas(srcdir=os.path.join(N1_outdir, processed_seuqences[monomer_sequence], 'DeepMSA2_a3m', 'finalMSAs'),
                                    trgdir=N1_monomer_deepmsa_outdir,
                                    srcname=processed_seuqences[monomer_sequence],
                                    trgname=monomer_id, rename_prefix=False)

            N1_monomer_outdir_img = os.path.join(N1_outdir_img, monomer_id)
            makedir_if_not_exists(N1_monomer_outdir_img)

            N2_monomer_outdir = os.path.join(N2_outdir, monomer_id)
            makedir_if_not_exists(N2_monomer_outdir)

            # make a copy
            N3_monomer_outdir = os.path.join(N3_outdir, monomer_id)

            if not os.path.exists(N3_monomer_outdir):
                os.system("cp -r " + os.path.join(N3_outdir, processed_seuqences[monomer_sequence]) + " " + N3_monomer_outdir)

    print("#################################################################################################")

    print("#################################################################################################")
    print("4. Start to generate complex alignments")

    N4_outdir = os.path.join(FLAGS.output_dir, 'N4_monomer_alignments_concatenation')
    makedir_if_not_exists(N4_outdir)
    
    is_homomers = len(processed_seuqences) == 1
    try:
        concat_methods = ['pdb_interact', 'species_interact', 'uniclust_oxmatch',
                           'string_interact', 'uniprot_distance', 'deepmsa2']
        if is_homomers:
            concat_methods = ['pdb_interact', 'species_interact', 'uniclust_oxmatch', 'deepmsa2'] 
        run_monomer_msas_concatenation_pipeline(
            multimer=','.join([chain_id for chain_id in chain_id_map]),
            run_methods=concat_methods,
            monomer_aln_dir=N1_outdir, monomer_model_dir=N3_outdir, outputdir=N4_outdir, params=params, is_homomers=is_homomers)
    except Exception as e:
        print(e)
        print("Program failed in step 5")

    print("#################################################################################################")

    print("#################################################################################################")

    print("5. Start to search complex templates based on monomer structures")

    N5_outdir = os.path.join(FLAGS.output_dir, 'N5_monomer_templates_concatenation')

    run_monomer_templates_concatenation_pipeline(multimers=[chain_id for chain_id in chain_id_map],
                                         monomer_aln_dir=N1_outdir,
                                         monomer_model_dir=N3_outdir,
                                         outdir=N5_outdir, params=params)

    print("#################################################################################################")

    print("#################################################################################################")

    print("6. Start to generate complex multimer structures")
    N6_outdir = os.path.join(FLAGS.output_dir, 'N6_multimer_structure_generation')

    makedir_if_not_exists(N6_outdir)
    run_methods = ['default_multimer']
    if not run_multimer_structure_generation_pipeline_v2(params=params,
                                                         fasta_path=FLAGS.fasta_path,
                                                         chain_id_map=chain_id_map,
                                                         aln_dir=N1_outdir,
                                                         complex_aln_dir=N4_outdir,
                                                         template_dir=N5_outdir,
                                                         monomer_model_dir=N3_outdir,
                                                         output_dir=N6_outdir,
                                                         run_methods=run_methods,
                                                         run_script=True,
                                                         run_deepmsa=False):
        print("Program failed in step 7")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
