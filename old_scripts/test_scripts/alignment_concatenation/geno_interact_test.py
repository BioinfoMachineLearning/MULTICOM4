import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from multicom_dev.tool import hhblits
from multicom_dev.tool import jackhmmer
from multicom_dev.common.util import is_dir, is_file, read_option_file, makedir_if_not_exists
from multicom_dev.complex.geno_dist import Geno_interact
from multicom_dev.alignment.alignment import *
from multicom_dev.complex.complex import write_concatenated_alignment

def run_hhblits(inparams):

    fasta, outdir, hhblits_binary, database = inparams

    hhblits_runner = hhblits.HHBlits(binary_path=hhblits_binary, databases=[database])

    return hhblits_runner.query(fasta, outdir)


def run_jackhmmer(inparams):

    fasta, outdir, jackhmmer_binary, database = inparams

    jackhmmer_runner = jackhmmer.Jackhmmer(binary_path=jackhmmer_binary, database_path=database)

    return jackhmmer_runner.query(fasta, outdir)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=is_file, required=True)
    parser.add_argument('--fasta1', type=is_file, required=True)
    parser.add_argument('--fasta2', type=is_file, required=True)
    parser.add_argument('--hhblits', type=is_file, required=True)
    parser.add_argument('--jackhmmer', type=is_file, required=True)
    parser.add_argument('--outdir', type=is_dir, required=True)

    args = parser.parse_args()

    params = read_option_file(args.option_file)

    # test hhblits
    outdir = args.outdir + '/hhblits'
    makedir_if_not_exists(outdir)

    process_list = []
    process_list.append([args.fasta1, outdir, args.hhblits, params['uniref_db_dir'] + '/' + params['uniref_db']])
    process_list.append([args.fasta2, outdir, args.hhblits, params['uniref_db_dir'] + '/' + params['uniref_db']])
    pool = Pool(processes=2)
    results = pool.map(run_hhblits, process_list)
    pool.close()
    pool.join()

    with open(results[0]['a3m']) as f:
        aln_1 = Alignment.from_file(f, format="a3m")

    with open(results[1]['a3m']) as f:
        aln_2 = Alignment.from_file(f, format="a3m")

    Geno_interact_runner = Geno_interact(params['uniprot_to_embl_table'], params['ena_genome_location_table'])
    pair_ids = Geno_interact_runner.get_interactions(aln_1, aln_2)

    target_header, sequences_full, sequences_monomer_1, sequences_monomer_2 = write_concatenated_alignment(pair_ids, aln_1, aln_2)

    # save the alignment files
    mon_alignment_file_1 = f"{outdir}/{aln_1.main_id}_monomer_1.a3m"
    with open(mon_alignment_file_1, "w") as of:
        write_a3m(sequences_monomer_1, of)

    mon_alignment_file_2 = f"{outdir}/{aln_2.main_id}_monomer_2.a3m"
    with open(mon_alignment_file_2, "w") as of:
        write_a3m(sequences_monomer_2, of)

    pair_ids.to_csv(f"{outdir}/{aln_1.main_id}_{aln_2.main_id}_interact.csv", index=False)
    print(pair_ids)

    complex_ailgnment_file = f"{outdir}/{aln_1.main_id}_{aln_2.main_id}.a3m"
    with open(complex_ailgnment_file, "w") as of:
        write_a3m(sequences_full, of)

    # test jackhmmer
    outdir = args.outdir + '/jackhmmer'
    makedir_if_not_exists(outdir)

    process_list = []
    process_list.append([args.fasta1, outdir, args.jackhmmer, params['uniref90_fasta']])
    process_list.append([args.fasta2, outdir, args.jackhmmer, params['uniref90_fasta']])
    pool = Pool(processes=2)
    results = pool.map(run_jackhmmer, process_list)
    pool.close()
    pool.join()

    with open(results[0]['sto']) as f:
        aln_1 = Alignment.from_file(f, format="stockholm")

    with open(results[1]['sto']) as f:
        aln_2 = Alignment.from_file(f, format="stockholm")

    pair_ids = Geno_interact_runner.get_interactions(aln_1, aln_2)

    target_header, sequences_full, sequences_monomer_1, sequences_monomer_2 =  write_concatenated_alignment(pair_ids, aln_1, aln_2)

    # save the alignment files
    mon_alignment_file_1 = f"{outdir}/{aln_1.main_id}_monomer_1.a3m"
    with open(mon_alignment_file_1, "w") as of:
        write_a3m(sequences_monomer_1, of)

    mon_alignment_file_2 = f"{outdir}/{aln_2.main_id}_monomer_2.a3m"
    with open(mon_alignment_file_2, "w") as of:
        write_a3m(sequences_monomer_2, of)

    pair_ids.to_csv(f"{outdir}/{aln_1.main_id}_{aln_2.main_id}_interact.csv", index=False)
    print(pair_ids)

    complex_ailgnment_file = f"{outdir}/{aln_1.main_id}_{aln_2.main_id}.a3m"
    with open(complex_ailgnment_file, "w") as of:
        write_a3m(sequences_full, of)
