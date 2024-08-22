import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle


class Alphafold_pkl_qa:

    def __init__(self):

        self.methods = ['ptm', 'iptm', 'plddt_avg', 'confidence', 'af3_ranking_score']

    def run(self, input_dir):
        ranking_pd = pd.DataFrame(columns=['model', 'ptm', 'iptm', 'plddt_avg', 'confidence', 'af3_ranking_score'])
        model_count = 0
        for pkl in os.listdir(input_dir):
            if pkl.find('.pkl') < 0 or pkl == 'features.pkl':
                continue
            ranking = {'model': pkl.replace('.pkl', '.pdb')}
            with open(os.path.join(input_dir, pkl), 'rb') as f:
                prediction_result = pickle.load(f)
                ranking['plddt_avg'] = np.mean(prediction_result['plddt'])
                ranking['ptm'] = float(prediction_result['ptm'])
                ranking['iptm'] = float(prediction_result['iptm'])
                ranking['confidence'] = float(prediction_result['ranking_confidence'])
                if 'ranking_score' in prediction_result:
                    ranking['af3_ranking_score'] = float(prediction_result['ranking_score'])
                else:
                    ranking['af3_ranking_score'] = -1
                ranking_pd = ranking_pd.append(pd.DataFrame(ranking, index=[model_count]))
                model_count += 1
        return ranking_pd.sort_values(by=['confidence'], ascending=False, ignore_index=True)
