import os
import pandas as pd

def get_model_scores(label_path):
    df = pd.read_csv(label_path)
    return {row['model']: row['tmscore'] for _, row in df.iterrows()}

def get_af_ranking_df(target_name, QA_dir):
    ranking_path = os.path.join(QA_dir, target_name.replace('-D1', ''), 'alphafold_ranking.csv')
    return pd.read_csv(ranking_path)

def extract_top5_models(ranking_df):
    """Returns top-5 models for each method."""
    method_models = {
        'af3': [],
        'af2': [],
    }

    for model in ranking_df['model']:
        base_name = model.replace('.pdb', '')
        if 'af3_jian' not in model and 'af3_pawan' not in model:
            if len(method_models['af2']) < 5:
                method_models['af2'].append(base_name)

    ranking_df = ranking_df.sort_values(by='af3_ranking_score', ascending=False).reset_index(drop=True)
    for model in ranking_df['model']:
        base_name = model.replace('.pdb', '')
        if 'af3_jian' in model or 'af3_pawan' in model:
            if len(method_models['af3']) < 5:
                method_models['af3'].append(base_name)

    return method_models

def best_score_among_top5(scores_dict, model_list):
    return max([scores_dict.get(m, -1) for m in model_list], default=-1)

if __name__ == "__main__":
    results_dir = "./target_result/GDT-TS/label"
    QA_dir = "./QAs/"

    # --- Section 1: AF2 vs AF3 (Best of Top-5) ---
    print("AF2 vs AF3 (Best of Top-5 models)")
    print("Target\tAF2\tAF3")

    for file in os.listdir(results_dir):
        if not file.endswith(".csv"):
            continue

        target = file.replace('.csv', '')
        label_path = os.path.join(results_dir, file)

        try:
            scores = get_model_scores(label_path)
            ranking_df = get_af_ranking_df(target, QA_dir)
            model_lists = extract_top5_models(ranking_df)

            af2_score = best_score_among_top5(scores, model_lists['af2'])
            af3_score = best_score_among_top5(scores, model_lists['af3'])

            print(f"{target}\t{af2_score}\t{af3_score}")

        except Exception as e:
            print(f"Error processing {target} in AF2 vs AF3: {e}")

    # Section 2: Alignment Comparison
    results_dir = "./target_result/GDT-TS/raw_result"
    print("\nAlignments: af3, colabfold, deepmsa, af2, dhr, esm-msa")
    print("Target\tColabFold\tDeepMSA_dMSA\tDeepMSA_qMSA\tDefault\tESM-MSA\tDHR")

    # Mapping of methods to their file names
    msa_methods = {
        # "af3": "af3.results",
        "colabfold": "colabfold_web.results",
        "deepmsa_dMSA": "deepmsa_dMSA.results",
        "deepmsa_qMSA": "deepmsa_qMSA.results",
        "default": "default.results",
        "esm_msa": "def_esm_msa.results",
        "dhr": "dhr.results"
    }

    # Collect scores per target per method
    all_scores = {}

    for method, filename in msa_methods.items():
        filepath = os.path.join(results_dir, filename)
        if not os.path.exists(filepath):
            print(f"Warning: {filename} not found.")
            continue
        with open(filepath) as f:
            for line in f:
                contents = line.strip().split()
                if len(contents) < 2:
                    continue
                target, top1_score = contents[0], max(contents[1:])
                if target not in all_scores:
                    all_scores[target] = {}
                all_scores[target][method] = top1_score

    # Print all scores in the expected column order
    for target in sorted(all_scores.keys()):
        scores = all_scores[target]
        print(f"{target}\t" +
              # f"{scores.get('af3', '-')}\t" +
              f"{scores.get('colabfold', '-')}\t" +
              f"{scores.get('deepmsa_dMSA', '-')}\t" +
              f"{scores.get('deepmsa_qMSA', '-')}\t" +
              f"{scores.get('default', '-')}\t" +
              f"{scores.get('esm_msa', '-')}\t" +
              f"{scores.get('dhr', '-')}")
