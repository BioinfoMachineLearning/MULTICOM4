import os, sys, argparse, time, copy, pathlib
from multiprocessing import Pool
from multicom4.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from multicom4.monomer_structure_refinement import iterative_refine_pipeline
from multicom4.common.protein import read_qa_txt_as_df, parse_fasta, make_chain_id_map, complete_result
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
    
    with open(FLAGS.fasta_path) as f:
        input_fasta_str = f.read()
    input_seqs, input_descs = parse_fasta(fasta_string=input_fasta_str)
    chain_id_map, chain_id_seq_map = make_chain_id_map(sequences=input_seqs, descriptions=input_descs)

    outdir = FLAGS.output_dir

    makedir_if_not_exists(outdir)


    N1_outdir = os.path.join(FLAGS.output_dir, 'N1_monomer_alignments_generation') 
    N3_outdir = os.path.join(FLAGS.output_dir, 'N3_monomer_structure_generation')

    N6_outdir = os.path.join(outdir, 'N6_multimer_structure_generation')
    
    default_workdir = os.path.join(N6_outdir, 'default_multimer')
    ranking_json_file = os.path.join(default_workdir, "ranking_debug.json")
    if not os.path.exists(ranking_json_file):
        raise Exception(f"Haven't generated default_multimer models!")

    run_methods = ["alphafold", "apollo", "bfactor"]

    N7_outdir = os.path.join(FLAGS.output_dir, 'N7_monomer_only_structure_evaluation')
    monomer_qas_res = {}
    processed_seuqences = {}
    for chain_id_idx, chain_id in enumerate(chain_id_map):
        monomer_id = chain_id
        monomer_sequence = chain_id_map[chain_id].sequence
        if monomer_sequence not in processed_seuqences:

            N1_monomer_outdir = os.path.join(N1_outdir, monomer_id)
            N7_monomer_outdir = os.path.join(N7_outdir, monomer_id)
            makedir_if_not_exists(N7_monomer_outdir)
            
            contact_map_file = os.path.join(N1_monomer_outdir, 'dncon4', f'{monomer_id}.dncon2.rr')
            if not os.path.exists(contact_map_file):
                raise Exception("The contact map file hasn't been generated!")

            dist_map_file = os.path.join(N1_monomer_outdir, 'deepdist', f'{monomer_id}.txt')
            if not os.path.exists(dist_map_file):
                raise Exception("The distance map file hasn't been generated!")

            result = run_monomer_evaluation_pipeline(params=params,
                                                    targetname=monomer_id,
                                                    fasta_file=os.path.join(FLAGS.output_dir, f"{monomer_id}.fasta"),
                                                    input_monomer_dir=os.path.join(N3_outdir, monomer_id),
                                                    input_multimer_dir="",
                                                    contact_map_file=contact_map_file,
                                                    dist_map_file=dist_map_file,
                                                    outputdir=N7_monomer_outdir, 
                                                    generate_final_models=False,
                                                    run_methods=run_methods)
            if result is None:
                raise RuntimeError(f"Program failed in step 7: monomer {monomer_id} model evaluation")
            monomer_qas_res[monomer_id] = result
            processed_seuqences[monomer_sequence] = monomer_id
        else:
            # make a copy
            N7_monomer_outdir = os.path.join(N7_outdir, monomer_id)
            os.system("cp -r " + os.path.join(N7_outdir, processed_seuqences[monomer_sequence]) + " " + N7_monomer_outdir)
            for msa in os.listdir(os.path.join(N7_monomer_outdir, 'msa')):
                os.system(f"sed -i 's/>{processed_seuqences[monomer_sequence]}/>{monomer_id}/g' " + os.path.join(N7_monomer_outdir, 'msa', msa))
            monomer_qas_res[monomer_id] = copy.deepcopy(monomer_qas_res[processed_seuqences[monomer_sequence]])


    bash_script_dir = os.path.join(N6_outdir, 'post_def_bash_scripts')
    if os.path.exists(params['slurm_script_template']):
        bash_script_dir = os.path.join(N6_outdir, 'post_def_slurm_scripts')
    os.makedirs(bash_script_dir, exist_ok=True)

    run_methods = ['folds_iter', 'folds_iter_not', 'folds_iter_o', 'folds_iter_o_not',
                   'folds_iter_esm', 'folds_iter_esm_not', 'folds_iter_esm_o', 'folds_iter_esm_o_not',
                   'def_mul_refine']

    for run_method in run_methods:

        if run_method == "def_mul_refine":
            
            for i in range(HETEROMULTIMER_CONFIG.predictors.def_mul_refine.number_of_input_models):

                cmd = f"python bin/multimer/heteromer_refine.py --option_file {FLAGS.option_file} " \
                f"--fasta_path {FLAGS.fasta_path} --output_dir {FLAGS.output_dir} " \
                f"--config_name {run_method} --idx {i}"

                if os.path.exists(params['slurm_script_template']):
                    bash_file = os.path.join(bash_script_dir, f"{run_method}_{i}.sh")
                    print(f"Generating bash file for {run_method}: {bash_file}")
                    jobname = f"{targetname}_{run_method}"
                    with open(bash_file, 'w') as fw:
                        for line in open(params['slurm_script_template']):
                            line = line.replace("JOBNAME", jobname)
                            fw.write(line)
                        fw.write(cmd)
                    #os.system(f"sbatch {bash_file}")
                else:
                    bash_file = os.path.join(bash_script_dir, f"{run_method}_{i}.sh")
                    print(bash_file)
                    with open(bash_file, 'w') as fw:
                        fw.write(cmd)

        else:

            cmd = f"python bin/multimer/homomer_foldseek.py --option_file {FLAGS.option_file} " \
                f"--fasta_path {FLAGS.fasta_path} --output_dir {FLAGS.output_dir} " \
                f"--config_name {run_method}"

            if os.path.exists(params['slurm_script_template']):
                bash_file = os.path.join(bash_script_dir, run_method + '.sh')
                print(f"Generating bash file for {run_method}: {bash_file}")
                jobname = f"{targetname}_{run_method}"
                with open(bash_file, 'w') as fw:
                    for line in open(params['slurm_script_template']):
                        line = line.replace("JOBNAME", jobname)
                        fw.write(line)
                    fw.write(cmd)
                os.system(f"sbatch {bash_file}")
            else:
                bash_file = os.path.join(bash_script_dir, run_method + '.sh')
                print(bash_file)
                with open(bash_file, 'w') as fw:
                    fw.write(cmd)
    

if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
