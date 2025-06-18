import os
import shutil
import subprocess
import argparse
import re
import pandas as pd

def run_tmscore(model_file, reference_file, tm_score_path, output_file):
    """Run TM-score and return the output."""
    if os.path.exists(output_file):
        print(f"Skipping TM-score for {model_file}, result already exists.")
        return

    # Including the original parameters -ter 1 -tmscore 6 for TM-score execution
    cmd = f"{tm_score_path} -ter 1 -tmscore 6 {model_file} {reference_file} > {output_file}"
    subprocess.run(cmd, shell=True, check=True)
    print(f"TM-score output saved to {output_file}")

def parse_tmscore_output(tm_result_file):
    """Parse TM-score output to extract the TM-score."""
    tm_score = None
    with open(tm_result_file, 'r') as f:
        for line in f:
            if line.startswith("TM-score") and "length of Structure_2" in line:
                tm_score = float(line.split()[1])
                break
    return tm_score

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

def evaluate_top_models(models_dir, target_id, temp_dir, pdb_filtered_file, workdir):
    """Evaluate the top models and return the TM-scores organized by model number."""
    scores_by_group = {}

    # af2 confidence score ranking
    alphafold_ranking_file = os.path.join(workdir, 'N7_multimer_structure_evaluation', 'alphafold_ranking.csv')
    af_ranking_df = pd.read_csv(alphafold_ranking_file)
    
    for i in range(len(af_ranking_df)):
        model_name = af_ranking_df.loc[i, 'model']
        if model_name.find('af3') >= 0:
            continue
        break
    
    tmscore_result_file = os.path.join(temp_dir, f"{model_name}_filtered_out")

    tm_score = parse_tmscore_output(tm_result_file)
    if tm_score is not None:
        group_name = 'af2'
        # Group ID is derived from the model name before the first underscore
        if model_name.find('af3_jian') == 0 or model_name.find('af3_pawan') == 0:
            group_name = 'af3'
        if group_name not in scores_by_group:
            scores_by_group[group_name] = []  # Initialize list for top_n models

        # Store the TM-score for the corresponding model number (i-1 index)
        scores_by_group[group_name] += [tm_score]

    print(f"Model: {model_name}, TM-score: {tm_score}")

    return scores_by_group

def evaluate_casp_models(pdb_filtered_dir, casp_models_dir, temp_dir, output_folder, ref_workdir_list, target_results_dir):
    """Evaluate CASP models for a list of targets."""
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if not os.path.exists(target_results_dir):
        os.makedirs(target_results_dir)

    # Clear target results directory
    for file in os.listdir(target_results_dir):
        os.remove(os.path.join(target_results_dir, file))

    # Read target list
    target_list = sorted(os.listdir(temp_dir))

    target_workdir_dict = {}
    for line in open(ref_workdir_list):
        target_id, workdir = line.rstrip('\n').split()
        target_workdir_dict[target_id] = target_workdir_dict

    all_results = {}

    for target_id in target_list:
        print(f"Evaluating Target {target_id}")
        temp_target_dir = os.path.join(temp_dir, target_id)
        if not os.path.exists(temp_target_dir):
            os.makedirs(temp_target_dir)

        pdb_filtered_file = os.path.join(pdb_filtered_dir, f"{target_id}_filtered.pdb")

        # Evaluate models for the target
        scores_by_group = evaluate_top_models(casp_models_dir, target_id, temp_target_dir, pdb_filtered_file, target_workdir_dict[target_id])

        # Save results for the target
        target_result_file = os.path.join(target_results_dir, f"{target_id}.results")
        with open(target_result_file, 'w') as result_file:
            #result_file.write("Group\tTM-score1\tTM-score2\tTM-score3\tTM-score4\tTM-score5\n")
            models, true_scores = [], []
            csv_file_path = os.path.join(output_folder, f"{target_id}.csv")
            for group in sorted(scores_by_group.keys()):
                scores = scores_by_group[group]
                result_file.write(f"{group:<20}" + "\t".join(f"{score:.4f}" if score is not None else "0.00" for score in scores) + "\n")

                group_result_file = os.path.join(target_results_dir, f"{group}.results")
                with open(group_result_file, 'a') as result_file2:
                    # Write target_id followed by scores
                    result_file2.write(f"{target_id:<20}" + " ".join(f"{score:.4f}" if score is not None else "0.00" for score in scores) + "\n")

                
                for i, score in enumerate(scores, start=1):
                    if score is None:
                        continue
                    if target_id.find('o') > 0:
                        true_target_id = target_id.replace('o', '')
                        models += [f"{true_target_id}{group}_{i}o"]
                    else:
                        models += [f"{target_id}{group}_{i}"]
                    true_scores += [float(score)]

            df = pd.DataFrame({'model': models, 'tmscore': true_scores})
            df.to_csv(csv_file_path)

        # Store results for summary
        all_results[target_id] = scores_by_group

def main():
    parser = argparse.ArgumentParser(description="Evaluate CASP models using TM-score.")
    parser.add_argument("--filtered_pdb_dir", help="Directory containing filtered PDB files.")
    parser.add_argument("--casp_models_dir", help="Directory containing CASP model predictions.")
    parser.add_argument("--temp_dir", help="Temporary directory for storing intermediate files.")
    parser.add_argument("--output_folder", help="Output directory for final ranking files.")
    parser.add_argument("--ref_workdir_list", help="File containing a list of target IDs to evaluate.")
    parser.add_argument("--target_results_dir", help="Directory for saving target-specific results.")

    args = parser.parse_args()

    evaluate_casp_models(
        args.filtered_pdb_dir,
        args.casp_models_dir,
        args.temp_dir,
        args.output_folder,
        args.ref_workdir_list,
        args.target_results_dir
    )

if __name__ == "__main__":
    main()
