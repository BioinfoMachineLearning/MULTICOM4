import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm

def correct_domain(domain_info_file, disorder_pred):
    contents = open(domain_info_file).readlines()

    domain_num = len(contents)
    print(f"Total {domain_num} domains")

    domain_disorder_info, domain_range_info = [],[]

    for line in contents:
        line = line.rstrip('\n')
        domain_range = line.split(':')[1]
        start, end = domain_range.split('-')

        start, end = int(start), int(end)

        disorder_region = disorder_pred[start-1:end]
        disorder_num = len([char for char in disorder_region if char == "T"])
        ratio = disorder_num / float(len(disorder_region))
        if ratio > 0.7:
            domain_disorder_info += ['Disorder']
        else:
            domain_disorder_info += ['Normal']

        domain_range_info += [f"{start}-{end}"]

    return domain_disorder_info, domain_range_info


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #parser.add_argument('--option_file', type=str, required=True)
    parser.add_argument('--fasta_path', type=str, required=True)
    parser.add_argument('--outputdir', type=str, required=True)
    parser.add_argument('--inpdb', type=str, required=False, default="")

    args = parser.parse_args()

    # params = read_option_file(args.option_file)

    os.makedirs(args.outputdir, exist_ok=True)

    # predict the disorder region
    diso_out = os.path.join(args.outputdir, 'disopred3')
    os.makedirs(diso_out, exist_ok=True)

    fasta_name = os.path.basename(args.fasta_path)

    diso_out_file = os.path.join(diso_out, fasta_name.replace('.fasta', '.diso'))
    if not os.path.exists(diso_out_file):
        os.system(f"cp {args.fasta_path} {diso_out}/{fasta_name}")
        print(f"perl /bmlfast/bml_casp16/tools/disopred/run_disopred.pl {diso_out}/{fasta_name} &> {diso_out}/run.log")
        os.system(f"perl /bmlfast/bml_casp16/tools/disopred/run_disopred.pl {diso_out}/{fasta_name} &> {diso_out}/run.log")

    sequence = open(args.fasta_path).readlines()[1].rstrip('\n')
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
    print(f"Predicting domain boundries using hhsearch")
    dom_hhsearch_out = os.path.join(args.outputdir, 'Dom_by_hhsearch')
    domain_info_file = os.path.join(dom_hhsearch_out, 'domain_info')
    if os.path.exists(domain_info_file):
        print(f"Found {domain_info_file}!!")
    else:
        hhsearchdb = "/bmlfast/bml_casp15/tools/casp14/Human_TS_CASP14_dependency/databases/hhsearch15db"
        print(f"perl /bmlfast/bml_casp15/tools/casp14/Human_TS_CASP14_dncon4_v2/scripts/P1_get_domains_by_hhsearch15.pl {args.fasta_path} {dom_hhsearch_out} {hhsearchdb} 40\n")
        os.system(f"perl /bmlfast/bml_casp15/tools/casp14/Human_TS_CASP14_dncon4_v2/scripts/P1_get_domains_by_hhsearch15.pl {args.fasta_path} {dom_hhsearch_out} {hhsearchdb} 40")
        os.system(f"cp {dom_hhsearch_out}/domain_info_thre40 {domain_info_file}")

    if os.path.exists(domain_info_file):
        domain_disorder_info, domain_range_info = correct_domain(domain_info_file, disorder_pred)
        domain_def = os.path.join(dom_hhsearch_out, 'domain_info_final')
        with open(domain_def, 'w') as fw:
            print("Corrected domain regions:")
            for i in range(len(domain_range_info)):
                print(f"domain {i}: {domain_range_info[i]} {domain_disorder_info[i]}")
                fw.write(f"domain {i}: {domain_range_info[i]} {domain_disorder_info[i]}\n")


    # 2. domain_parser, requires the predicted full-length model as input
    dom_parse_out = os.path.join(args.outputdir, 'Dom_by_domain_parser')
    domain_info_file = os.path.join(dom_parse_out, 'domain_info')
    if len(args.inpdb) > 0 and os.path.exists(args.inpdb):
        print(f"Found input pdb, using domain parser to predict domain boundries")
        os.makedirs(dom_parse_out, exist_ok=True)
        if os.path.exists(domain_info_file):
            print(f"Found {domain_info_file}!!")
        else:
            print(f"perl /bmlfast/bml_casp15/tools/casp14/Human_TS_CASP14_dncon4_v2/scripts/P1_get_domains_by_DomainParser_inpdb.pl {args.fasta_path} {args.inpdb} {dom_parse_out} ")
            os.system(f"perl /bmlfast/bml_casp15/tools/casp14/Human_TS_CASP14_dncon4_v2/scripts/P1_get_domains_by_DomainParser_inpdb.pl {args.fasta_path} {args.inpdb} {dom_parse_out} ")

    if os.path.exists(domain_info_file):
        domain_disorder_info, domain_range_info = correct_domain(domain_info_file, disorder_pred)
        domain_def = os.path.join(dom_parse_out, 'domain_info_final')
        with open(domain_def, 'w') as fw:
            print("Corrected domain regions:")
            for i in range(len(domain_range_info)):
                print(f"domain {i}: {domain_range_info[i]} {domain_disorder_info[i]}")
                fw.write(f"domain {i}: {domain_range_info[i]} {domain_disorder_info[i]}\n")


    # 3. Unidoc structure-based
    unidoc_parse_out = os.path.join(args.outputdir, 'Dom_by_UniDoc')
    domain_info_file = os.path.join(unidoc_parse_out, 'domain_info')
    if len(args.inpdb) > 0 and os.path.exists(args.inpdb):
        print(f"Found input pdb, using UniDoc to predict domain boundries")
        os.makedirs(unidoc_parse_out, exist_ok=True)
        if os.path.exists(domain_info_file):
            print(f"Found {domain_info_file}!!")
        else:
            os.chdir("/bmlfast/bml_casp16/tools/UniDoc")
            print(f"python Run_UniDoc_from_scratch_structure.py -i {args.inpdb} -c A > {unidoc_parse_out}/unidoc.out")
            os.system(f"python Run_UniDoc_from_scratch_structure.py -i {args.inpdb} -c A > {unidoc_parse_out}/unidoc.out")

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
            
