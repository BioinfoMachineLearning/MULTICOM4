import os, sys
import pandas as pd
import numpy as np
from multicom4.common.pipeline import run_monomer_msa_pipeline

def run_hhsearch_dom(disorder_pred, N1_outdir, fasta_path, params):
    # 1. hhsearch
    print(f"Predicting domain boundries using hhsearch")
    dom_hhsearch_out = os.path.join(N1_outdir, 'Dom_by_hhsearch')
    domain_info_file = os.path.join(dom_hhsearch_out, 'domain_info')
    if os.path.exists(domain_info_file):
        print(f"Found {domain_info_file}!!")
    else:
        cmd = f"perl {params['dom_hhsearch_script']} {fasta_path} {dom_hhsearch_out} {params['hhsearch15db']} 40"
        print(cmd)
        os.system(cmd)
        os.system(f"cp {dom_hhsearch_out}/domain_info_thre40 {domain_info_file}")

    if os.path.exists(domain_info_file):
        domain_disorder_info, domain_range_info = correct_domain(domain_info_file, disorder_pred)
        domain_def = os.path.join(dom_hhsearch_out, 'domain_info_final')
        with open(domain_def, 'w') as fw:
            print("Corrected domain regions:")
            for i in range(len(domain_range_info)):
                print(f"domain {i}: {domain_range_info[i]} {domain_disorder_info[i]}")
                fw.write(f"domain {i}: {domain_range_info[i]} {domain_disorder_info[i]}\n")


def run_dom_parser(disorder_pred, N1_outdir, fasta_path, ranked_0_pdb, params):
    dom_parse_out = os.path.join(N1_outdir, 'Dom_by_domain_parser')
    domain_info_file = os.path.join(dom_parse_out, 'domain_info')
    print(f"Found input pdb, using domain parser to predict domain boundries")
    if len(sequence) < 3000:     
        os.makedirs(dom_parse_out, exist_ok=True)
        if os.path.exists(domain_info_file):
            print(f"Found {domain_info_file}!!")
        else:
            cmd = f"perl {params['dom_parser_script']} {fasta_path} {ranked_0_pdb} {dom_parse_out}"
            print(cmd)
            os.system(cmd)
    else:
        print(f"Domain parser cannot deal with sequence longer than 3000, skip!")

    if os.path.exists(domain_info_file):
        domain_disorder_info, domain_range_info = correct_domain(domain_info_file, disorder_pred)
        domain_def = os.path.join(dom_parse_out, 'domain_info_final')
        with open(domain_def, 'w') as fw:
            print("Corrected domain regions:")
            for i in range(len(domain_range_info)):
                print(f"domain {i}: {domain_range_info[i]} {domain_disorder_info[i]}")
                fw.write(f"domain {i}: {domain_range_info[i]} {domain_disorder_info[i]}\n")


def run_unidoc(disorder_pred, N1_outdir, fasta_path, ranked_0_pdb, params):
    unidoc_parse_out = os.path.join(N1_outdir, 'Dom_by_UniDoc')
    domain_info_file = os.path.join(unidoc_parse_out, 'domain_info')
    print(f"Found input pdb, using UniDoc to predict domain boundries")
    os.makedirs(unidoc_parse_out, exist_ok=True)
    if os.path.exists(domain_info_file):
        print(f"Found {domain_info_file}!!")
    else:
        os.chdir(params['unidoc_program_path'])
        cmd = f"python {params['unidoc_program_script']} -i {ranked_0_pdb} -c A > {unidoc_parse_out}/unidoc.out"
        print(cmd)
        os.system(cmd)

        domain_range_out = open(os.path.join(unidoc_parse_out, 'unidoc.out')).readlines()[0]

        domain_ranges = domain_range_out.rstrip('\n').split('/')

        sorted_domain_ranges = sorted(domain_ranges, key=lambda x: int(x.split('~')[0]))
        print(sorted_domain_ranges)

        with open(domain_info_file, 'w') as fw:
            for i, domain_str in enumerate(sorted_domain_ranges):
                domain_str = domain_str.replace('~', '-')
                fw.write(f"domain {i}: {domain_str}\n")

    if os.path.exists(domain_info_file):
        domain_disorder_info, domain_range_info = correct_domain(domain_info_file, disorder_pred)
        domain_def = os.path.join(unidoc_parse_out, 'domain_info_final')
        with open(domain_def, 'w') as fw:
            print("Corrected domain regions:")
            for i in range(len(domain_range_info)):
                print(f"domain {i}: {domain_range_info[i]} {domain_disorder_info[i]}")
                fw.write(f"domain {i}: {domain_range_info[i]} {domain_disorder_info[i]}\n")


def generate_domain_fastas(fasta_path, domain_info, output_dir):

    os.makedirs(outputdir, exist_ok=True)

    sequence = open(fasta_path).readlines()[1].rstrip('\n')
    for line in open(domain_info):
        line = line.rstrip('\n')
        # domain 1: 140-317 Normal
        # domain 1: 140-317
        domain_name, domain_range = line.split(':')
        domain_name = domain_name.replace(' ', '')
        domain_range = domain_range.split()[0]

        domain_sequence = ""
        starts, ends = [], []

        for domain_range_str in domain_range.split(','):
            start, end = domain_range_str.split('-')
            starts += [int(start)-1]
            ends += [int(end)]

        for start, end in zip(starts, ends):
            domain_sequence += sequence[start:end]

        with open(os.path.join(outputdir, domain_name + '.fasta'), 'w') as fw:
            fw.write(f">{domain_name}\n{domain_sequence}")


def generate_domain_alignments(params, domain_dir, output_dir):

    makedir_if_not_exists(output_dir)

    for fasta_path in os.listdir(domain_dir):

        fasta_path = os.path.join(domain_dir, fasta_path)

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

        outdir = os.path.join(output_dir, targetname)

        makedir_if_not_exists(outdir)

        N1_outdir = os.path.join(outdir, 'monomer_alignments')
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
