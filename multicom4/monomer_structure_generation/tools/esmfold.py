import os, sys, argparse, time
import torch
import esm

def get_bfactors(infile):
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
    return np.array(bfactors)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    # parser.add_argument('--num_recycle', type=int, required=False, default=4)
    parser.add_argument('--chunk_size', type=int, required=False, default=None)

    args = parser.parse_args()

    sequences = []
    for line in open(args.fasta):
        if line[0] == ">":
            continue
        sequences += [line.rstrip('\n')]
    
    sequence = sequences[0]
    if len(sequences) > 1:
        sequence = ':'.join(sequences)
    
    try:
        for num_recycle in range(4, 104, 2):

            model = esm.pretrained.esmfold_v1()
            model = model.eval().cuda()
            model.set_chunk_size(args.chunk_size)

            with torch.no_grad():
                output = model.infer_pdb(sequence, num_recycles=num_recycle)

            outpdb = os.path.join(args.outdir, f"{num_recycle}.pdb")
            with open(outpdb, "w") as f:
                f.write(output)

            # import biotite.structure.io as bsio
            # struct = bsio.load_structure("result.pdb", extra_fields=["b_factor"])
            # print(struct.b_factor.mean())  # this will be the pLDDT
            # # 88.3
    except Exception as e:
        print(e)

    if len(sequences) == 1: 
        ranking_confidences = {}
        ranked_order = []
        for pdbfile in os.listdir(args.outdir):
            if pdbfile.find('.pdb') < 0:
                continue
            
            pdbname = pdbfile.replace('.pdb', '')

            plddts = get_bfactors(os.path.join(args.outdir, pdbfile))

            result_output_path = os.path.join(args.outdir, f'result_{pdbname}.pkl')
            with open(result_output_path, 'wb') as f:
                pickle.dump({'plddt': plddts}, f, protocol=4)

            ranked_order.append(pdbname)
            ranking_confidences[pdbname] = np.mean(plddts)

        ranked_order = [model_name for model_name, confidence in sorted(ranking_confidences.items(), key=lambda x: x[1], reverse=True)]

        # Write out relaxed PDBs in rank order.
        for idx, model_name in enumerate(ranked_order):
            ranked_output_path = os.path.join(args.outdir, f'ranked_{idx}.pdb')
            os.system(f"cp {os.path.join(args.outdir, model_name + '.pdb')} {ranked_output_path}")

        ranking_output_path = os.path.join(args.outdir, 'ranking_debug.json')
        with open(ranking_output_path, 'w') as f:
            f.write(json.dumps({'plddts': ranking_confidences, 'order': ranked_order}, indent=4))


