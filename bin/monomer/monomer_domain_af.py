import os, sys, argparse, time, copy, pathlib
from multiprocessing import Pool
from multicom4.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from multicom4.monomer_domain_combination.domaim import run_hhsearch_dom, run_dom_parser, run_unidoc
import pandas as pd
from absl import flags
from absl import app

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to monomer fasta')
flags.DEFINE_string('output_dir', None, 'Output directory')
flags.DEFINE_string('domain_info', None, 'Output directory')
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

    fasta_name = os.path.basename(FLAGS.fasta_path)
    diso_out_file = os.path.join(N1_outdir, fasta_name.replace('.fasta', '.diso'))
    if not os.path.exists(diso_out_file):
        raise Exception(f"Disorder prediction hasn't been generated!")

    disorder_pred = []
    for line in open(diso_out_file):
        if line.startswith('#'):
            continue
        tmps = line.split()
        if tmps[2] == "*":
            disorder_pred += ["T"]
        else:
            disorder_pred += ["N"]

    print("Disorder prediction:")
    print(sequence)
    print(''.join(disorder_pred))

    # 1. hhsearch
    run_hhsearch_dom(disorder_pred=disorder_pred, 
                     N1_outdir=N1_outdir,
                     fasta_path=FLAGS.fasta_path,
                     params=params)

    N3_outdir = os.path.join(outdir, 'N3_monomer_structure_generation')    
    default_workdir = os.path.join(N3_outdir, 'default')
    ranking_json_file = os.path.join(default_workdir, "ranking_debug.json")
    if not os.path.exists(ranking_json_file):
        raise Exception(f"Haven't generated default models!")

    ranked_0_pdb = os.path.join(N3_outdir, 'default', 'ranked_0.pdb')
    # 2. domain_parser, requires the predicted full-length model as input
    run_dom_parser(disorder_pred=disorder_pred, 
                   N1_outdir=N1_outdir,
                   fasta_path=FLAGS.fasta_path,
                   ranked_0_pdb=ranked_0_pdb,
                   params=params)


    # 3. Unidoc structure-based
    run_unidoc(disorder_pred=disorder_pred, 
               N1_outdir=N1_outdir,
               fasta_path=FLAGS.fasta_path,
               ranked_0_pdb=ranked_0_pdb,
               params=params)

    




























    bash_script_dir = os.path.join(N3_outdir, 'post_def_bash_scripts')
    if os.path.exists(params['slurm_script_template']):
        bash_script_dir = os.path.join(N3_outdir, 'post_def_slurm_scripts')
    os.makedirs(bash_script_dir, exist_ok=True)

    run_methods = ['default_tmsearch', 'def_esm_msa', 'foldseek_refine', 'foldseek_refine_esm', 'foldseek_refine_esm_h']

    for run_method in run_methods:
        cmd = f"python bin/monomer/monomer_single_predictor.py --option_file {FLAGS.option_file} " \
              f"--fasta_path {FLAGS.fasta_path} --output_dir {FLAGS.output_dir} " \
              f"--config_name {run_method}"
        print(f"Generating bash scripts for {run_method}")

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