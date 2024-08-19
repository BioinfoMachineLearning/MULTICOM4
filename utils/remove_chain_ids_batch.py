import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    parser.add_argument('--chain_ids', type=str, required=True)
    args = parser.parse_args()

    remove_ids = args.chain_ids.split(',')

    for inpdb in os.listdir(args.indir):
        print(f"Processing {inpdb}")
        chain_contents = {}
        for line in open(args.indir + '/' + inpdb):
            if line.startswith('ATOM'):
                chain_id = line[21]
                if chain_id not in chain_contents:
                    chain_contents[chain_id] = [line]
                else:
                    chain_contents[chain_id] += [line]

        with open(args.outdir + '/' + inpdb, 'w') as fw:
            for chain_id in chain_contents:
                if chain_id in remove_ids:
                    continue
                for line in chain_contents[chain_id]:
                    fw.write(line)
                fw.write("TER\n")

