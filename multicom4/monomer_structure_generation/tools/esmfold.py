import os, sys, argparse, time
import torch
import esm

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=str, required=True)
    parser.add_argument('--outpdb', type=str, required=True)
    parser.add_argument('--num_recycle', type=int, required=False, default=4)
    parser.add_argument('--chunk_size', type=int, required=False, default=None)

    args = parser.parse_args()

    sequence = open(args.fasta).readlines()[1].rstrip('\n')

    try:
        model = esm.pretrained.esmfold_v1()
        model = model.eval().cuda()
        model.set_chunk_size(args.chunk_size)

        with torch.no_grad():
            output = model.infer_pdb(sequence, num_recycles=args.num_recycle)

        with open(args.outpdb, "w") as f:
            f.write(output)

        # import biotite.structure.io as bsio
        # struct = bsio.load_structure("result.pdb", extra_fields=["b_factor"])
        # print(struct.b_factor.mean())  # this will be the pLDDT
        # # 88.3
    except Exception as e:
        print(e)


