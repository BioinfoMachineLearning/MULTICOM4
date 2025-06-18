import os
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Evaluate CASP models using TM-score.")
    parser.add_argument("--rankingdir", required=True, help="Directory containing filtered PDB files.")
    parser.add_argument("--labeldir", required=True, help="Directory containing filtered PDB files.")
    args = parser.parse_args()

    for targetname in sorted(os.listdir(args.rankingdir)):
        if targetname[0] != 'T' and targetname[0] != 'H':
            continue

        alphafold_ranking_file = os.path.join(args.rankingdir, targetname, 'alphafold_ranking.csv')
        af_ranking_df = pd.read_csv(alphafold_ranking_file)

        af_ranking_df = af_ranking_df.sort_values(by=['confidence'], ascending=False)
        af_ranking_df.reset_index(inplace=True, drop=True)
        top1_model_af_confidence = af_ranking_df.loc[0, 'model']

        af_ranking_df = af_ranking_df.sort_values(by=['af3_ranking_score'], ascending=False)
        af_ranking_df.reset_index(inplace=True, drop=True)
        top1_model_af3_ranking = af_ranking_df.loc[0, 'model']

        pairwise_ranking_file = os.path.join(args.rankingdir, targetname, 'multieva.csv')
        pairwise_ranking_df = pd.read_csv(pairwise_ranking_file)
        pairwise_ranking_df = pairwise_ranking_df.sort_values(by=['MMalign score'], ascending=False)
        pairwise_ranking_df.reset_index(inplace=True, drop=True)
        top1_model_pairwise_ranking = pairwise_ranking_df.loc[0, 'Name']

        gate_ranking_file = os.path.join(args.rankingdir, targetname, 'gate.csv')
        gate_ranking_df = pd.read_csv(gate_ranking_file)
        gate_ranking_df = gate_ranking_df.sort_values(by=['score'], ascending=False)
        gate_ranking_df.reset_index(inplace=True, drop=True)
        top1_model_gate_ranking = gate_ranking_df.loc[0, 'model']

        gate_score_dict = {}
        for i in range(len(gate_ranking_df)):
            model = gate_ranking_df.loc[i, 'model']
            gate_score = gate_ranking_df.loc[i, 'score']
            gate_score_dict[model] = gate_score

        af3_gate_avg_dict = {'model': [], 'avg_score': []}
        for i in range(len(af_ranking_df)):
            model = af_ranking_df.loc[i, 'model']
            ranking_score = af_ranking_df.loc[i, 'af3_ranking_score']
            if ranking_score < 0:
                continue
            af3_gate_avg_dict['model'] += [model]
            af3_gate_avg_dict['avg_score'] += [(ranking_score + gate_score_dict[model]) / 2]
        af3_gate_avg_df = pd.DataFrame(af3_gate_avg_dict)
        af3_gate_avg_df = af3_gate_avg_df.sort_values(by=['avg_score'], ascending=False)
        af3_gate_avg_df.reset_index(inplace=True, drop=True)
        # top1_model_avg_ranking = af3_gate_avg_df.loc[0, 'model']

        af2_gate_avg_dict = {'model': [], 'avg_score': []}
        for i in range(len(af_ranking_df)):
            model = af_ranking_df.loc[i, 'model']
            ranking_score = af_ranking_df.loc[i, 'confidence']
            af2_gate_avg_dict['model'] += [model]
            af2_gate_avg_dict['avg_score'] += [(ranking_score + gate_score_dict[model]) / 2]
        af2_gate_avg_df = pd.DataFrame(af2_gate_avg_dict)
        af2_gate_avg_df = af2_gate_avg_df.sort_values(by=['avg_score'], ascending=False)
        af2_gate_avg_df.reset_index(inplace=True, drop=True)
        top1_model_avg_ranking = af2_gate_avg_df.loc[0, 'model']

        target_labels_df = pd.read_csv(os.path.join(args.labeldir, targetname + '.csv'))
        # print(targetname)
        score_dict = {target_labels_df.loc[i, 'model'] + '.pdb' : target_labels_df.loc[i, 'tmscore'] for i in range(len(target_labels_df))}

        #top1_scores = [targetname]
        #top1_scores += [str(score_dict[top1_model_af_confidence])]
        #top1_scores += [str(score_dict[top1_model_af3_ranking])]
        #top1_scores += [str(score_dict[top1_model_pairwise_ranking])]
        #top1_scores += [str(score_dict[top1_model_gate_ranking])]
        top1_scores = [str(score_dict[top1_model_avg_ranking])]

        print('\t'.join(top1_scores))

if __name__ == "__main__":
    main()




