import os
import shutil
import subprocess
import argparse
import re
import pandas as pd
import json

def parse_tmscore_output(tm_result_file):
    # print(tm_result_file)
    with open(tm_result_file) as f:
        data = json.load(f)
        return data

def find_model_file(models_dir, target_id, model_number):
    """Find the model file by checking multiple naming conventions."""
    target_dir = os.path.join(models_dir, target_id)
    if not os.path.exists(target_dir):
        print(f"Target directory {target_dir} does not exist.")
        return []

    # Allow for variations in target names (e.g., with or without "o" or other suffixes)
    if target_id.endswith("o"):
        base_target_id = target_id[:-1]  # Remove the "o" for pattern matching
    else:
        base_target_id = target_id

    # Regular expression to match patterns like T1201TS272_2 or T1201TS272_2o
    pattern = re.compile(rf"{base_target_id}TS\d+_{model_number}o?_filtered.pdb")

    matching_files = [
        os.path.join(target_dir, filename)
        for filename in os.listdir(target_dir)
        if pattern.match(filename)
    ]
    return matching_files

def evaluate_top_models(models_dir, target_id, temp_dir, pdb_filtered_file, top_n = 10):
    """Evaluate the top models and return the TM-scores organized by model number."""
    keys = ['ICS', 'IPS', 'QSbest', 'DockQ_wave', 'LDDT', 'TMscore_os']
    all_scores_by_group = {'ICS': {}, 'IPS': {}, 'QSbest': {}, 'DockQ_wave': {}, 'LDDT': {}, 'TMscore_os': {}}
    for i in range(1, top_n + 1):
        # Find matching model files
        # matching_files = find_model_file(models_dir, target_id, i)
        # if not matching_files:
        #     print(f"Model file for {target_id} with model number {i} does not exist, skipping.")
        #     continue
        print(sorted(os.listdir(os.path.join(models_dir, target_id))))
        matching_files = [os.path.join(models_dir, target_id, model_file) for model_file in sorted(os.listdir(os.path.join(models_dir, target_id))) if int(model_file.split('_')[-2][0]) == i-1 ]

        for model_file in matching_files:
            model_name = os.path.basename(model_file)
            model_filtered = os.path.join(temp_dir, f"{model_name}_filtered")
            shutil.copyfile(model_file, model_filtered)
            tm_result_file = f"{model_filtered}_out"

            if not os.path.exists(tm_result_file):
                raise Exception(f"Cannot find {tm_result_file}")
               
            # Parse TM-score
            result_dict = parse_tmscore_output(tm_result_file)
            # group_id = "TS" + model_name.split("_")[0].split("TS")[1]
            # group_id = "_".join(model_name.split("_")[:2])
            group_id = "_".join(model_name.split("_")[:-2])
            print(group_id)
            for key in all_scores_by_group:
                if group_id not in all_scores_by_group[key]:
                    all_scores_by_group[key][group_id] = [None] * top_n

            try:
                all_scores_by_group['ICS'][group_id][i-1] = float(result_dict['ics_trimmed'])
            except Exception:
                all_scores_by_group['ICS'][group_id][i-1] = 0.0

            try:
                all_scores_by_group['IPS'][group_id][i-1] = float(result_dict['ips_trimmed'])
            except Exception:
                all_scores_by_group['IPS'][group_id][i-1] = 0.0

            try:
                all_scores_by_group['QSbest'][group_id][i-1] = float(result_dict['qs_best'])
            except Exception:
                all_scores_by_group['QSbest'][group_id][i-1] = 0.0

            try:
                all_scores_by_group['DockQ_wave'][group_id][i-1] = float(result_dict['dockq_wave'])
            except Exception:
                all_scores_by_group['DockQ_wave'][group_id][i-1] = 0.0

            try:
                all_scores_by_group['LDDT'][group_id][i-1] = float(result_dict['lddt'])
            except Exception:
                all_scores_by_group['LDDT'][group_id][i-1] = 0.0

            try:
                all_scores_by_group['TMscore_os'][group_id][i-1] = float(result_dict['tm_score'])
            except Exception:
                all_scores_by_group['TMscore_os'][group_id][i-1] = 0.0

    return all_scores_by_group

def evaluate_casp_models(pdb_filtered_dir, casp_models_dir, temp_dir, target_list_file, target_results_dir):
    """Evaluate CASP models for a list of targets."""

    if not os.path.exists(target_results_dir):
        os.makedirs(target_results_dir)

    # Read target list
    with open(target_list_file, 'r') as f:
        target_list = [line.strip() for line in f if line.strip()]

    all_results = {}

    for target_id in target_list:
        print(f"Evaluating Target {target_id}")
        temp_target_dir = os.path.join(temp_dir, target_id)
        if not os.path.exists(temp_target_dir):
            os.makedirs(temp_target_dir)

        pdb_filtered_file = os.path.join(pdb_filtered_dir, f"{target_id}_filtered.pdb")

        # Evaluate models for the target
        all_scores_by_group = evaluate_top_models(
            casp_models_dir, target_id, temp_target_dir, pdb_filtered_file, 10
        )

        # print(all_scores_by_group['TMscore_os'])

        for score_name in all_scores_by_group:
            result_dir = os.path.join(target_results_dir, score_name)
            os.makedirs(result_dir, exist_ok=True)

            raw_dir = os.path.join(result_dir, 'raw_result')
            os.makedirs(raw_dir, exist_ok=True)

            label_dir = os.path.join(result_dir, 'label')
            os.makedirs(label_dir, exist_ok=True)

            # Save results for the target
            target_result_file = os.path.join(raw_dir, f"{target_id}.results")
            with open(target_result_file, 'w') as result_file:
                #result_file.write("Group\tTM-score1\tTM-score2\tTM-score3\tTM-score4\tTM-score5\n")
                models, true_scores = [], []
                csv_file_path = os.path.join(label_dir, f"{target_id}.csv")
                for group in sorted(all_scores_by_group[score_name].keys()):
                    scores = all_scores_by_group[score_name][group]
                    # print(scores)
                    result_file.write(f"{group:<20}" + "\t".join(f"{score:.4f}" if score is not None else "0.00" for score in scores) + "\n")

                    group_result_file = os.path.join(raw_dir, f"{group}.results")
                    with open(group_result_file, 'a') as result_file2:
                        # Write target_id followed by scores
                        result_file2.write(f"{target_id:<20}" + " ".join(f"{score:.4f}" if score is not None else "0.00" for score in scores) + "\n")

                    
                    for i, score in enumerate(scores):
                        if score is None:
                            continue
                        models += [f"{group}_{i}"]
                        true_scores += [float(score)]

                df = pd.DataFrame({'model': models, 'tmscore': true_scores})
                df.to_csv(csv_file_path)
        

def main():
    parser = argparse.ArgumentParser(description="Evaluate CASP models using TM-score.")
    parser.add_argument("--filtered_pdb_dir", help="Directory containing filtered PDB files.")
    parser.add_argument("--casp_models_dir", help="Directory containing CASP model predictions.")
    parser.add_argument("--temp_dir", help="Temporary directory for storing intermediate files.")
    parser.add_argument("--target_list_file", help="File containing a list of target IDs to evaluate.")
    parser.add_argument("--target_results_dir", help="Directory for saving target-specific results.")

    args = parser.parse_args()

    evaluate_casp_models(
        args.filtered_pdb_dir,
        args.casp_models_dir,
        args.temp_dir,
        args.target_list_file,
        args.target_results_dir,
    )

if __name__ == "__main__":
    main()
