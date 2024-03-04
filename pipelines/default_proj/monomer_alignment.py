import os, sys, argparse, time, copy, pathlib
from multiprocessing import Pool
from multicom4.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from multicom4.common.protein import read_qa_txt_as_df, complete_result
from multicom4.common.pipeline import run_monomer_msa_pipeline, run_monomer_template_search_pipeline
import pandas as pd
from absl import flags
from absl import app

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fastadir', None, 'Path to multimer fastas')
flags.DEFINE_string('output_dir', None, 'Output directory')
FLAGS = flags.FLAGS

def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    check_file(FLAGS.option_file)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    params['colabfold_search_program'] = ""
    params['colabfold_split_msas_program'] = ""
    params['colabfold_databases'] = ""
    params['mmseq_program'] = ""

    params['deepmsa2_path'] = ""
    params['JGIclust_database'] = ""
    params['metaclust_database'] = ""
    
    params['DHR_program_path'] = ""
    params['DHR_database_path'] = ""

    for fasta_path in os.listdir(FLAGS.fastadir):

        fasta_path = FLAGS.fastadir + '/' + fasta_path

        # targetname = pathlib.Path(FLAGS.fasta_path).stem
        targetname = None
        sequence = None
        for line in open(fasta_path):
            line = line.rstrip('\n').strip()
            if line.startswith('>'):
                targetname = line[1:].split()[0]
                # if targetname_in_fasta != targetname:
                #     print("Warning: fasta file name doesn't match with fasta content!")
            else:
                sequence = line

        outdir = FLAGS.output_dir + '/' + targetname

        makedir_if_not_exists(outdir)

        N1_outdir = outdir + '/N1_monomer_alignments_generation'
        makedir_if_not_exists(N1_outdir)

        print("#################################################################################################")
        print(f"1. Start to generate alignments for monomers")

        uniref30 = params['uniref_db']
        uniclust30 = ''
        uniref90_fasta = params['uniref90_fasta']

        smallbfd = ""  # params['smallbfd_database']
        bfd = params['bfd_database']
        mgnify = params['mgnify_database']

        hhblits_binary = params['hhblits_program']
        hhfilter_binary = params['hhfilter_program']
        jackhmmer_binary = params['jackhmmer_program']

        colabfold_search_binary = '' #params['colabfold_search_program']
        colabfold_split_msas_binary = '' #params['colabfold_split_msas_program']
        colabfold_databases = '' #params['colabfold_databases']
        mmseq_binary = '' #params['mmseq_program']

        deepmsa2_path = '' #params['deepmsa2_path']
        JGIclust_database_path = '' #params['JGIclust_database']
        metaclust_database_path = '' #params['metaclust_database']
        
        dhr_program_path = '' #params['DHR_program_path']
        dhr_database_path = '' #params['DHR_database_path']


        result = run_monomer_msa_pipeline(fasta=fasta_path, outdir=N1_outdir, params=params, only_monomer=True)

        if result is None:
            raise RuntimeError('The monomer alignment generation has failed!')

        # N1_outdir_img = outdir + '/N1_monomer_alignments_generation_img'
        # makedir_if_not_exists(N1_outdir_img)
        # img_msa = run_monomer_msa_pipeline_img(params=params, fasta=fasta_path, outdir=N1_outdir_img)

        # print("#################################################################################################")
        # print("2. Start to generate template for monomer")

        # N2_outdir = outdir + '/N2_monomer_template_search'

        # makedir_if_not_exists(N2_outdir)

        # template_file = run_monomer_template_search_pipeline(params=params, targetname=targetname, sequence=sequence,
        #                                                     a3m=f"{N1_outdir}/{targetname}_uniref90.sto", outdir=N2_outdir)

        # if template_file is None:
        #     raise RuntimeError("Program failed in step 2: monomer template search")

        # print("The generation for monomer template has finished!")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fastadir',
        'output_dir'
    ])
    app.run(main)
