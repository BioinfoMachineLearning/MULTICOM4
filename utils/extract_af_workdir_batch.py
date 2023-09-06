import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
from multicom_dev.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file, is_file, is_dir


def extract_scores_from_pkl(inpkl, outpkl):
    with open(inpkl, 'rb') as f:
        prediction_result = pickle.load(f)  
        prediction_result_new = {}
        prediction_result_new['iptm'] = prediction_result['iptm']
        prediction_result_new['plddt'] = prediction_result['plddt']
        prediction_result_new['ptm'] = prediction_result['ptm']
        prediction_result_new['ranking_confidence'] = prediction_result['ranking_confidence']
        prediction_result_new['max_predicted_aligned_error'] = prediction_result['max_predicted_aligned_error']
        prediction_result_new['predicted_aligned_error'] = prediction_result['predicted_aligned_error']
        with open(outpkl, 'wb') as f:
            pickle.dump(prediction_result_new, f, protocol=4)
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=is_file, required=True)
    args = parser.parse_args()

    for target in os.listdir(args.indir):
        af_workdir = args.indir + '/' + target + '/N6_quaternary_structure_generation'
        if not os.path.exists(af_workdir):
            raise Exception(f"Cannot find the {af_workdir}!")
        
        af_scores_workdir = args.indir + '/' + target + '/N6_only_scores'
        if os.path.exists(af_scores_workdir):
            os.system(f"rm -rf {af_scores_workdir}")
        makedir_if_not_exists(af_scores_workdir)

        for method in os.listdir(af_workdir):
            srcdir = af_workdir + '/' + method
            trgdir = af_scores_workdir + '/' + method
            makedir_if_not_exists(trgdir)

            for af_file in os.listdir(srcdir):
                if af_file.find('.pkl') < 0 or af_file == 'features.pkl':
                    os.system(f"cp -r {srcdir}/{af_file} {trgdir}")
                else:
                    extract_scores_from_pkl(srcdir + '/' + af_file, 
                                            trgdir + '/' + af_file)
