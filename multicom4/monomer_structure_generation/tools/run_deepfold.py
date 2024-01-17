import os, sys, argparse, time
import numpy as np
import json

def get_avg_factor(infile):
    bfactors = []
    read = set()
    for line in open(infile):
        if line.startswith('ATOM'):
            resn = line[17:20]
            chain = line[21]
            resi = line[22:26]
            icode = line[26]
            bfactor = line[60:66]
            r_uid = (resn, chain, resi, icode)
            if r_uid not in read:
                read.add(r_uid)
                #print(read)
            else:
                continue
            bfactors += [float(bfactor)]
    return np.mean(np.array(bfactors))

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--pickle_path', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    # parser.add_argument('--deepfold_dir', type=str, required=True)
    args = parser.parse_args()

    # os.chdir(args.deepfold_dir)
    os.chdir('/bmlfast/bml_casp16/tools/deepfold')
    cmds = []
    models = ['model1', 'model2', 'model3', 'model4', 'model5']
    for model in models:
        relaxed_model = f'{args.outdir}/prot_00000/relaxed_{model}.pdb'
        if os.path.exists(relaxed_model):
                continue
        cmd = f"python run_from_pkl.py " \
              f"--pickle_paths {args.pickle_path} " \
              f"--model_names {model} " \
              f"--model_paths params/{model}.npz " \
              f"--output_dir {args.outdir} "
        cmds += [cmd]

    for cmd in cmds:
        try:
            print(cmd)
            os.system(cmd)
        except Exception as e:
            print(e)

    ranking_confidences = {}

    for modelnum in range(1, 6):
        relaxed_model = f'{args.outdir}/prot_00000/relaxed_model{modelnum}.pdb'
        pklfile = f'{args.outdir}/prot_00000/result_model{modelnum}.pkl'

        modelname = f"model_{modelnum}_pred_0"
        os.system(f"cp {relaxed_model} {args.outdir}/relaxed_{modelname}.pdb")
        os.system(f"cp {pklfile} {args.outdir}/result_{modelname}.pkl")

        ranking_confidences[modelname] =  get_avg_factor(relaxed_model)
    
    # Rank by model confidence.
    ranked_order = [model_name for model_name, confidence in sorted(ranking_confidences.items(), key=lambda x: x[1], reverse=True)]

    # Write out relaxed PDBs in rank order.
    for idx, model_name in enumerate(ranked_order):
        ranked_output_path = os.path.join(args.outdir, f'ranked_{idx}.pdb')
        os.system(f"cp {args.outdir}/relaxed_{model_name}.pdb {ranked_output_path}")

    ranking_output_path = os.path.join(args.outdir, 'ranking_debug_deepfold.json')
    with open(ranking_output_path, 'w') as f:
        f.write(json.dumps({'plddts': ranking_confidences, 'order': ranked_order}, indent=4))
