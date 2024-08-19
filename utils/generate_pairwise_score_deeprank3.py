import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--targetname', type=str, required=True)
    parser.add_argument('--infile', type=str, required=True)
    parser.add_argument('--outfile', type=str, required=True)

    args = parser.parse_args()

    data_dict = {'model': [], 'score': []}
    df = pd.read_csv(args.infile, index_col=[0])
    for model_idx, model in enumerate(df.columns):
        tmscores = [df[model][i] for i in range(len(df[model])) if i != model_idx]
        data_dict['model'] += [model]
        data_dict['score'] += [np.mean(np.array(tmscores))]
    
    pairwise_df = pd.DataFrame(data_dict)
    pairwise_df = pairwise_df.sort_values(by='score', ascending=False)
    pairwise_df.reset_index(inplace=True, drop=True)

    contents = ['PFRMAT QA']
    contents += [f"TARGET {args.targetname}"]
    contents += ['MODEL 1']
    contents += ['QMODE 1']
    
    for i in range(len(pairwise_df)):
        model = pairwise_df.loc[i, 'model']
        score = pairwise_df.loc[i, 'score']
        contents += [f"{model} {score}"]

    contents += ['END']

    with open(args.outfile, 'w') as fw:
        fw.write('\n'.join(contents))

