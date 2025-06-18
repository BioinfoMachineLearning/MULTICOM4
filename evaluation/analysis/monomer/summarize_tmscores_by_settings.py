import os
import pandas as pd

def get_model_scores(label_path):
    df = pd.read_csv(label_path)
    return {row['model']: row['tmscore'] for _, row in df.iterrows()}

def get_af_ranking_df(target_name, QA_dir):
    af_ranking_path = os.path.join(QA_dir, target_name.replace('-D1', ''), 'alphafold_ranking.csv')
    return pd.read_csv(af_ranking_path)

def safe_score(scores_dict, model_name):
    return scores_dict.get(model_name, -1)

if __name__ == "__main__":
    results_dir = "./target_result/GDT-TS/label"
    QA_dir = "./QAs/"

    # Section 1: AF2 vs AF3
    print("AF2 vs AF3")
    print("Target\tAF2\tAF3")

    for file in os.listdir(results_dir):
        target = file.replace('.csv', '')
        label_path = os.path.join(results_dir, file)
        scores = get_model_scores(label_path)
        ranking_df = get_af_ranking_df(target, QA_dir)
        af2_model = ""
        for model in ranking_df['model']:
            base_name = model.replace('.pdb', '')
            if 'af3_jian' not in model and 'af3_pawan' not in model:
                if not af2_model:
                    af2_model = base_name
                    break

        ranking_df = ranking_df.sort_values(by='af3_ranking_score', ascending=False).reset_index(drop=True)

        af3_model = ranking_df.loc[0, 'model'].replace('.pdb', '')
        # print(f"{target}\t{af2_model}\t{af3_model}")
        print(f"{target}\t{safe_score(scores, af2_model)}\t{safe_score(scores, af3_model)}")

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
                target, top1_score = contents[0], contents[1]
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
