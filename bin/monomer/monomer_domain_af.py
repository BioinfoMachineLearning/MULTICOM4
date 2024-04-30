import os, sys, argparse, time, copy, pathlib
from multiprocessing import Pool
from multicom4.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from multicom4.monomer_domain_combination.domain import *
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
    targetname = os.path.basename(FLAGS.fasta_path).replace('.fasta', '')

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
    N1_outdir_img = os.path.join(outdir, 'N1_monomer_alignments_generation_img')
    N2_outdir = os.path.join(outdir, 'N2_monomer_template_search')

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

    domain_info_dict = {}

    # 1. hhsearch
    domain_info_dict['dom_hhsearch'] = run_hhsearch_dom(disorder_pred=disorder_pred, 
                                                        N1_outdir=N1_outdir,
                                                        fasta_path=FLAGS.fasta_path,
                                                        params=params)
    
    N3_outdir = os.path.join(outdir, 'N3_monomer_structure_generation')    
    default_workdir = os.path.join(N3_outdir, 'default')
    ranking_json_file = os.path.join(default_workdir, "ranking_debug.json")
    if not os.path.exists(ranking_json_file):
        raise Exception(f"Haven't generated default models!")

    alphafold_def_a3m = os.path.join(N3_outdir, 'default', 'msas', 'monomer_final.a3m')

    ranked_0_pdb = os.path.join(N3_outdir, 'default', 'ranked_0.pdb')
    # 2. domain_parser, requires the predicted full-length model as input

    if len(sequence) >= 3000:
        print(f"Domain parser cannot deal with sequence longer than 3000, skip!")
    else:
        dom_parser_domain_file = run_dom_parser(disorder_pred=disorder_pred, 
                                                N1_outdir=N1_outdir,
                                                fasta_path=FLAGS.fasta_path,
                                                ranked_0_pdb=ranked_0_pdb,
                                                params=params)
        
        if len(dom_parser_domain_file) > 0:
            domain_info_dict['dom_parser'] = dom_parser_domain_file 

    # 3. Unidoc structure-based
    domain_info_dict['dom_unidoc'] = run_unidoc(disorder_pred=disorder_pred, 
                                                N1_outdir=N1_outdir,
                                                fasta_path=FLAGS.fasta_path,
                                                ranked_0_pdb=ranked_0_pdb,
                                                params=params)

    if FLAGS.domain_info is not None:

        domain_info_dict['dom_manual'] = FLAGS.domain_info

    for run_method in domain_info_dict:

        method_outdir = os.path.join(N3_outdir, run_method)
        os.makedirs(method_outdir, exist_ok=True)

        fastadir = os.path.join(method_outdir, 'domain_fastas')
        os.makedirs(fastadir, exist_ok=True)

        generate_domain_fastas(fasta_path=FLAGS.fasta_path, 
                               domain_info=domain_info_dict[run_method], 
                               output_dir=fastadir)
        
        alndir = os.path.join(method_outdir, 'domain_alns')
        os.makedirs(alndir, exist_ok=True)

        generate_domain_alignments(params=params, 
                                   domain_fastas_dir=fastadir,
                                   output_dir=alndir)
        
        paired_a3m = combine_domain_msas(sequence=sequence, 
                                        domain_info=domain_info_dict[run_method], 
                                        domaindir=alndir)
        
        combine_a3m = os.path.join(alndir, 'final.a3m')
        with open(combine_a3m, 'w') as fw:
            fw.write(open(alphafold_def_a3m).read())
            fw.write(open(paired_a3m).read())
    

    if not run_monomer_structure_generation_pipeline_v2(params=params,
                                                        targetname=targetname,
                                                        fasta_path=FLAGS.fasta_path,
                                                        alndir=N1_outdir, 
                                                        img_alndir=N1_outdir_img,
                                                        templatedir=N2_outdir, 
                                                        outdir=N3_outdir,
                                                        run_methods=list(domain_info_dict.keys()),
                                                        run_script=os.path.exists(params['slurm_script_template'])):

        print("Program failed in step 3: monomer structure generation")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
