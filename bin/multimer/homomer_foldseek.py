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
    run_multimer_structure_generation_homo_pipeline_v2, \
    run_multimer_structure_generation_pipeline_foldseek, run_multimer_structure_generation_pipeline_foldseek_old, \
    run_multimer_evaluation_pipeline, run_monomer_msa_pipeline_img, \
    foldseek_iterative_monomer_input, copy_same_sequence_msas

from absl import flags
from absl import app
import copy
import pandas as pd

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to multimer fasta')
flags.DEFINE_string('output_dir', None, 'Output directory')
flags.DEFINE_string('config_name', None, 'config name')
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
    input_seqs, input_descs = parse_fasta(input_fasta_str)
    chain_id_map, chain_id_seq_map = make_chain_id_map(sequences=input_seqs,
                                                       descriptions=input_descs)

    print("#################################################################################################")

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

    if FLAGS.config_name == "def_mul_refine":
        N6_outdir = os.path.join(FLAGS.output_dir, 'N6_multimer_structure_generation')
        if not run_multimer_structure_generation_homo_pipeline_v2(params=params,
                                                                fasta_path=FLAGS.fasta_path,
                                                                chain_id_map=chain_id_map,
                                                                aln_dir=N1_outdir,
                                                                complex_aln_dir=N4_outdir,
                                                                template_dir=N5_outdir,
                                                                monomer_model_dir=N3_outdir,
                                                                output_dir=N6_outdir,
                                                                run_methods=[FLAGS.config_name],
                                                                run_script=True,
                                                                run_deepmsa=False):
            print("Program failed in step 6")
            
        print("Multimer structure generation has been finished!")
    
    else:
        
        print("Multimer structure generation has been finished!")

        print("#################################################################################################")

        print("7. Start to evaluate monomer models")

        N7_outdir = os.path.join(FLAGS.output_dir, 'N7_monomer_only_structure_evaluation')

        print("#################################################################################################")

        print("8. Start to run multimer iterative generation pipeline using top-ranked monomer models")

        qa_result_dir = N7_outdir

        iterative_prepare_dir = os.path.join(qa_result_dir, 'iter_prepare')
        makedir_if_not_exists(iterative_prepare_dir)

        pipeline_inputs = []
        for i in range(2):
            monomer_pdb_dirs = {}
            monomer_alphafold_a3ms = {}
            pdb_name = None

            first_monomer_id = ""
            for chain_id in chain_id_map:
                first_monomer_id = chain_id
                ranking_file = os.path.join(N7_outdir, first_monomer_id, 'pairwise_ranking_monomer.csv')
                monomer_ranking = pd.read_csv(ranking_file)
                pdb_name = monomer_ranking.loc[i, 'model']
                break

            current_work_dir = os.path.join(iterative_prepare_dir, str(i))
            makedir_if_not_exists(current_work_dir)

            for chain_id in chain_id_map:
                monomer_id = chain_id
                chain_pdb_dir = os.path.join(current_work_dir, monomer_id)
                makedir_if_not_exists(chain_pdb_dir)

                first_monomer_id_dir = os.path.join(qa_result_dir, first_monomer_id)
                os.system("cp " + os.path.join(first_monomer_id_dir, 'pdb', pdb_name) + " " + os.path.join(chain_pdb_dir, pdb_name))

                new_contents = []
                for idx, line in enumerate(open(os.path.join(qa_result_dir, first_monomer_id, 'msa', pdb_name.replace('.pdb', '.a3m'))).readlines()):
                    if idx == 0:
                        new_contents += [f">{monomer_id}\n"]
                    else:
                        new_contents += [line]

                open(os.path.join(chain_pdb_dir, pdb_name.replace('.pdb', '.a3m')), 'w').writelines(new_contents)
                monomer_pdb_dirs[chain_id] = os.path.join(chain_pdb_dir, pdb_name)
                monomer_alphafold_a3ms[chain_id] = os.path.join(chain_pdb_dir, pdb_name.replace('.pdb', '.a3m'))

            print(monomer_alphafold_a3ms)
            pipeline_inputs += [foldseek_iterative_monomer_input(monomer_pdb_dirs=monomer_pdb_dirs,
                                                                monomer_alphafold_a3ms=monomer_alphafold_a3ms)]

        monomer_template_stos = []
        for chain_id in chain_id_map:
            monomer = chain_id
            monomer_template_sto = os.path.join(N1_outdir, monomer, f"{monomer}_uniref90.sto")
            if not os.path.exists(monomer_template_sto):
                raise Exception(f"Cannot find template stos for {monomer}: {monomer_template_sto}")
            monomer_template_stos += [monomer_template_sto]

        if FLAGS.config_name.find('_o') > 0:
            if not run_multimer_structure_generation_pipeline_foldseek_old(params=params, fasta_path=FLAGS.fasta_path,
                                                                            chain_id_map=chain_id_map, config_name=FLAGS.config_name,
                                                                            pipeline_inputs=pipeline_inputs, outdir=N6_outdir,
                                                                            is_homomers=True, monomer_template_stos=monomer_template_stos):
                print(f"Program failed in step 6 {FLAGS.config_name}")
        else:
            if not run_multimer_structure_generation_pipeline_foldseek(params=params, fasta_path=FLAGS.fasta_path,
                                                                        chain_id_map=chain_id_map, config_name=FLAGS.config_name,
                                                                        pipeline_inputs=[pipeline_inputs[0]],
                                                                        outdir=N6_outdir,
                                                                        is_homomers=True, monomer_template_stos=monomer_template_stos):
                print("Program failed in step 6 foldseek_iter")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
