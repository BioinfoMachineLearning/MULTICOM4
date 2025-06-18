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
    
    for i in range(len(af_ranking_df)):
        model_name = af_ranking_df.loc[i, 'model']
        if model_name.find('af3') >= 0:
            continue
        break

    tmscore_result_file = os.path.join(temp_dir, f"{model_name}_filtered.pdb_filtered_out")

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
