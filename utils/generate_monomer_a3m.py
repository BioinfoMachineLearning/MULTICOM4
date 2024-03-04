import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from multicom4.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
#from multicom4.complex_alignment_generation.pipeline.concatenate import Concatenate
from multicom4.tool import hhblits
from multicom4.tool import jackhmmer
from multicom4.monomer_alignment_generation.pipeline import Monomer_alignment_generation_pipeline
from multicom4.monomer_templates_concatenation import parsers

MAX_MSA_DEPTH = 10

def write_a3m(msa, outfile):
    with open(outfile, 'w') as fw:
        for desc, sequence in zip(msa.descriptions, msa.sequences):
            fw.write(f">{desc}\n{sequence}\n")


def generate_a3ms_for_single_seq(inparams):

    input_fasta_path, outdir, params = inparams

    uniref30_database_path = params['uniref_db']
    bfd_database_path = params['bfd_database']
    mgnify_database_path = params['mgnify_database']
    uniref90_database_path = params['uniref90_fasta']
    
    hhblits_binary = params['hhblits_program']
    jackhmmer_binary = params['jackhmmer_program']

    jackhmmer_uniref90_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary,
                database_path=uniref90_database_path,
                get_tblout=False)
    
    uniref90_a3m_depth = 0
    uniref90_a3m = os.path.join(outdir, 'uniref90.a3m')

    if os.path.exists(uniref90_a3m) and len(open(uniref90_a3m).readlines()) <= 1:
        os.system(f"rm {uniref90_a3m}")

    if not os.path.exists(uniref90_a3m):
        uniref_sto = os.path.join(outdir, 'uniref90.sto')
        if os.path.exists(uniref_sto) and len(open(uniref_sto).readlines()) <= 5:
            os.system(f"rm {uniref_sto}")

        if not os.path.exists(uniref_sto):
            result = jackhmmer_uniref90_runner.query(input_fasta_path, uniref_sto)
        with open(uniref_sto) as f:
            uniref90_msa = parsers.parse_stockholm(f.read())
            uniref90_a3m_depth = len(uniref90_msa)
            uniref90_msa = uniref90_msa.truncate(max_seqs=10000)
            write_a3m(uniref90_msa, uniref90_a3m) 
        os.system(f"rm {uniref_sto}")  
    else:
        with open(uniref90_a3m) as f:
            msa = parsers.parse_a3m(f.read())
            uniref90_a3m_depth = len(msa) - 1

    if uniref90_a3m_depth > MAX_MSA_DEPTH:
        os.system(f"touch {outdir}/STOP")
        return

    jackhmmer_mgnify_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary,
                database_path=mgnify_database_path,
                get_tblout=False)
    
    mgnify_a3m_depth = 0
    mgnify_a3m = os.path.join(outdir, 'mgnify.a3m')
    if os.path.exists(mgnify_a3m) and len(open(mgnify_a3m).readlines()) <= 1:
        os.system(f"rm {mgnify_a3m}")

    if not os.path.exists(mgnify_a3m):
        mgnify_sto = os.path.join(outdir, 'mgnify.sto')
        if os.path.exists(mgnify_sto) and len(open(mgnify_sto).readlines()) <= 5:
            os.system(f"rm {mgnify_sto}")

        if not os.path.exists(mgnify_sto):
            result = jackhmmer_mgnify_runner.query(input_fasta_path, mgnify_sto) 
        with open(mgnify_sto) as f:
            mgnify_msa = parsers.parse_stockholm(f.read())
            mgnify_a3m_depth = len(mgnify_msa)
            mgnify_msa = mgnify_msa.truncate(max_seqs=10000)
            write_a3m(mgnify_msa, mgnify_a3m)   
        os.system(f"rm {mgnify_sto}") 
    else:
        with open(mgnify_a3m) as f:
            msa = parsers.parse_a3m(f.read())
            mgnify_a3m_depth = len(msa) - 1

    if uniref90_a3m_depth + mgnify_a3m_depth> MAX_MSA_DEPTH:
        os.system(f"touch {outdir}/STOP")
        return

    uniref30_bfd_msa_runner = hhblits.HHBlits(binary_path=hhblits_binary,
                                              databases=[bfd_database_path, uniref30_database_path])

    uniref30_bfd_a3m = os.path.join(outdir, 'uniref30_bfd.a3m')
    if not os.path.exists(uniref30_bfd_a3m):
        result = uniref30_bfd_msa_runner.query(input_fasta_path, uniref30_bfd_a3m)

    uniref30_bfd_a3m_depth = 0
    with open(uniref30_bfd_a3m) as f:
        bfd_msa = parsers.parse_a3m(f.read())
        uniref30_bfd_a3m_depth = len(bfd_msa) - 1

    if uniref30_bfd_a3m_depth + uniref90_a3m_depth + mgnify_a3m_depth > MAX_MSA_DEPTH:
        os.system(f"touch {outdir}/STOP")
        return

    os.system(f"touch {outdir}/DONE")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=is_file, required=True)
    parser.add_argument('--fastadir', type=is_dir, required=True)
    parser.add_argument('--output_dir', type=str, required=True)

    args = parser.parse_args()

    params = read_option_file(args.option_file)

    makedir_if_not_exists(args.output_dir)

    check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    process_list = []
    for fastafile in os.listdir(args.fastadir):
        targetname = fastafile.replace('.fasta', '')
        outdir = os.path.join(args.output_dir, targetname)
        makedir_if_not_exists(outdir)
        if os.path.exists(outdir + '/DONE') or os.path.exists(outdir + '/STOP'):
            print(f"{fastafile} has been processed!")
            continue
        os.system(f"cp {args.fastadir}/{fastafile} {outdir}")
        process_list.append([f"{outdir}/{fastafile}", outdir, params])

    print(f"Total {len(process_list)} monomers to be processed")
    pool = Pool(processes=10)
    results = pool.map(generate_a3ms_for_single_seq, process_list)
    pool.close()
    pool.join()
