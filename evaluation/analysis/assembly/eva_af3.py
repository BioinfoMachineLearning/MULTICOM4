import os
import shutil
import subprocess
import argparse
import re
import pandas as pd
import pickle
import json

def parse_tmscore_output(tm_result_file):
    # print(tm_result_file)
    with open(tm_result_file) as f:
        data = json.load(f)
        return data

def check_same_contents(file1, file2):
    """Check if two PDB files have the same content."""
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        return f1.read() == f2.read()

def evaluate_top_models(models_dir, target_id, temp_dir, pdb_filtered_file, workdir):
    """Evaluate the top models and return the TM-scores organized by model number."""

    # af2 confidence score ranking
    alphafold_ranking_file = os.path.join(workdir, 'N7_multimer_structure_evaluation', 'alphafold_ranking.csv')
    af_ranking_df = pd.read_csv(alphafold_ranking_file)
    
    top1_af3_model = ""
    if target_id == "H0208":
        top1_af3_model = 'af3_pawan_0.pdb'
    if target_id == "H0215" or target_id == "T0206o":
        top1_af3_model = 'af3_jian_0.pdb'

    if top1_af3_model == "":
        if 'af3_ranking_score' in af_ranking_df.columns:
            best_score = -1
            for i in range(len(af_ranking_df)):
                ranking_score = af_ranking_df.loc[i, 'af3_ranking_score']
                if ranking_score > best_score:
                    best_score = ranking_score
                    top1_af3_model = af_ranking_df.loc[i, 'model']
        else:
            af3_jian_workdir = os.path.join(workdir, 'N6_multimer_structure_generation', 'af3_jian')
            af3_pawan_workdir = os.path.join(workdir, 'N6_multimer_structure_generation', 'af3_pawan')
            best_score = float('-inf')
            for af3_workdir in [af3_jian_workdir, af3_pawan_workdir]:
                for pkl in os.listdir(af3_workdir):
                    if pkl.find('.pkl') < 0 or pkl == 'features.pkl':
                        continue
                    pkl_file = os.path.join(af3_workdir, pkl)
                    with open(pkl_file, 'rb') as f:
                        prediction_result = pickle.load(f)
                        #print(prediction_result)
                        af3_ranking_score = float(prediction_result['ranking_score'])
                        if af3_ranking_score > best_score:
                            best_score = af3_ranking_score
                            top1_af3_model = pkl_file.replace('.pkl', '.pdb')

            for af3_pdb in os.listdir(models_dir):
                if af3_pdb.find('af3_jian') < 0 and af3_pdb.find('af3_pawan') < 0:
                    continue
                if check_same_contents(top1_af3_mode, os.path.join(models_dir, af3_pdb)):
                    top1_af3_model = af3_pdb

    model_name = top1_af3_model
   
    tmscore_result_file = os.path.join(temp_dir, f"{model_name}_filtered.pdb_filtered_out")

    # print(tmscore_result_file)

    result_dict = parse_tmscore_output(tmscore_result_file)
    
    return float(result_dict['ics_trimmed']), float(result_dict['ips_trimmed']), float(result_dict['qs_best']), float(result_dict['dockq_wave']), float(result_dict['lddt']), float(result_dict['tm_score'])

def evaluate_casp_models(pdb_filtered_dir, casp_models_dir, temp_dir, ref_workdir_list):
    """Evaluate CASP models for a list of targets."""

    # Read target list
    target_list = sorted(os.listdir(temp_dir))

    target_workdir_dict = {}
    for line in open(ref_workdir_list):
        workdir, target_id  = line.rstrip('\n').split()
        target_workdir_dict[target_id] = workdir

    all_results = {}

    for target_id in target_list:
        # print(f"Evaluating Target {target_id}")
        temp_target_dir = os.path.join(temp_dir, target_id)
        if not os.path.exists(temp_target_dir):
            os.makedirs(temp_target_dir)

        pdb_filtered_file = os.path.join(pdb_filtered_dir, f"{target_id}_filtered.pdb")

        # Evaluate models for the target
        ics, ips, qs_best, dockq_wave, lddt, tmscore = evaluate_top_models(casp_models_dir, target_id, temp_target_dir, pdb_filtered_file, target_workdir_dict[target_id])

        print(f"{target_id}\t{ics}\t{ips}\t{qs_best}\t{dockq_wave}\t{lddt}\t{tmscore}")


def main():
    parser = argparse.ArgumentParser(description="Evaluate CASP models using TM-score.")
    parser.add_argument("--filtered_pdb_dir", help="Directory containing filtered PDB files.")
    parser.add_argument("--casp_models_dir", help="Directory containing CASP model predictions.")
    parser.add_argument("--temp_dir", help="Temporary directory for storing intermediate files.")
    parser.add_argument("--ref_workdir_list", help="File containing a list of target IDs to evaluate.")

    args = parser.parse_args()

    evaluate_casp_models(
        args.filtered_pdb_dir,
        args.casp_models_dir,
        args.temp_dir,
        args.ref_workdir_list
    )

if __name__ == "__main__":
    main()
