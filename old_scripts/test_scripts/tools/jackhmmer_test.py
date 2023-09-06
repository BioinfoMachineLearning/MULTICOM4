import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from multicom_dev.tool import hhblits
from multicom_dev.tool import jackhmmer
from multicom_dev.common.util import is_dir, is_file, read_option_file


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=is_file, required=True)
    parser.add_argument('--fasta', type=is_file, required=True)
    parser.add_argument('--jackhmmer', type=is_file, required=True)
    parser.add_argument('--outdir', type=is_dir, required=True)

    args = parser.parse_args()

    params = read_option_file(args.option_file)

    jackhmmer_runner = jackhmmer.Jackhmmer(
        binary_path=args.jackhmmer,
        database_path=params['uniref90_fasta'])

    jackhmmer_result = jackhmmer_runner.query(args.fasta, args.outdir)

    print(jackhmmer_result)

