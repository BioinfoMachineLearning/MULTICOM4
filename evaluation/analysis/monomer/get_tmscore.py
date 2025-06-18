import os
import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Evaluate CASP models using TM-score.")
    parser.add_argument("--rankingdir", required=True, help="Directory containing filtered PDB files.")
    parser.add_argument("--labeldir", required=True, help="Directory containing filtered PDB files.")
    args = parser.parse_args()

    for targetname in sorted(os.listdir(args.rankingdir)):
        if targetname[0] != 'T' and targetname[0] != 'H':
            continue
        
        # print(targetname)
        alphafold_ranking_file = os.path.join(args.rankingdir, targetname, 'alphafold_ranking.csv')
        af_ranking_df = pd.read_csv(alphafold_ranking_file)

        top1_models = []

        af_ranking_df = af_ranking_df.sort_values(by=['plddt_avg'], ascending=False)
        af_ranking_df.reset_index(inplace=True, drop=True)
        top1_model_af_confidence = af_ranking_df.loc[0, 'model']
        top1_models += [top1_model_af_confidence]

        pairwise_ranking_file = os.path.join(args.rankingdir, targetname, 'pairwise_ranking_monomer.csv')
        pairwise_ranking_df = pd.read_csv(pairwise_ranking_file)
        pairwise_ranking_df = pairwise_ranking_df.sort_values(by=['score'], ascending=False)
        pairwise_ranking_df.reset_index(inplace=True, drop=True)
        top1_model_pairwise_ranking = pairwise_ranking_df.loc[0, 'model']
        top1_models += [top1_model_pairwise_ranking]

        gate_ranking_file = os.path.join(args.rankingdir, targetname, 'gate_af_summary.csv')
        gate_ranking_df = pd.read_csv(gate_ranking_file)

        for method in ["gate", "enqa", "gcpnet_ema", "deeprank3_cluster", "deeprank3_singleqa", "deeprank3_singleqa_lite"]:
            gate_ranking_df = gate_ranking_df.sort_values(by=method, ascending=False)
            gate_ranking_df.reset_index(inplace=True, drop=True)
            top1_model_gate_ranking = gate_ranking_df.loc[0, 'model']
            top1_models += [top1_model_gate_ranking]


        target_labels_df = pd.read_csv(os.path.join(args.labeldir, targetname + '-D1.csv'))
        # print(targetname)
        score_dict = {target_labels_df.loc[i, 'model'] + '.pdb' : target_labels_df.loc[i, 'tmscore'] for i in range(len(target_labels_df))}
        max_scores = str(np.max(target_labels_df['tmscore']))

        top1_scores = [targetname]

        for top1_model in top1_models:
            if top1_model not in score_dict:
                top1_model += ".pdb"
            top1_scores += [str(score_dict[top1_model])]
            
        print('\t'.join(top1_scores + [max_scores]))

if __name__ == "__main__":
    main()




