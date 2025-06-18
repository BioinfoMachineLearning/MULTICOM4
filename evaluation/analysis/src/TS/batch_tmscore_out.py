import os, sys, argparse, time
from multiprocessing import Pool

def run_command(inparams):
    mmalign_program, pdb1, pdb2, outfile = inparams
    cmd = f"{mmalign_program} {pdb1} {pdb2} > {outfile} "
    print(cmd)
    os.system(cmd)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str, required=True)
    parser.add_argument('--nativedir', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    parser.add_argument('--mmalign_program', type=str, required=True)

    args = parser.parse_args()

    process_list = []

    for native_pdb in os.listdir(args.nativedir):

        targetname = native_pdb.replace('.pdb', '')

        outdir = args.outdir + '/' + targetname

        modeldir = os.path.join(args.indir, targetname)
        if not os.path.exists(modeldir):
            continue
        
        os.makedirs(outdir, exist_ok=True)

        for model in os.listdir(modeldir):
            os.system(f"cp {args.indir}/{targetname}/{model} {outdir}")
            process_list.append([args.mmalign_program, args.nativedir + '/' + native_pdb, modeldir + '/' + model, outdir + '/' + model + '_out'])

    pool = Pool(processes=150)
    results = pool.map(run_command, process_list)
    pool.close()
    pool.join()
