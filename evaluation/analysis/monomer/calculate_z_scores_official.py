import os
import pandas as pd

targets = os.listdir('./TSdomains')

targets_filt = [target.replace('.pdb' , '') for target in targets]

# Function to calculate z-scores
def calculate_z_scores(df, column):
    mean = df[column].mean()
    std = df[column].std()
    if std == 0:
        return df[column] * 0  # Avoid division by zero
    return (df[column] - mean) / std

# Function to extract group name from model identifier
def extract_group_name(model_id):
    # Assuming the model_id format is 'H1202TS051_1'
    parts = model_id.split('TS')
    if len(parts) > 1:
        return 'TS' + parts[1].split('_')[0]
    return None

# Function to process and generate per-target Z-score tables (rows as targets, columns as groups)
def generate_per_target_z_score_tables(base_path, metrics):
    for metric in metrics:
        
        metric_results = []
        
        for file in os.listdir(base_path):  # Assuming CSV files in each folder
            if not file.endswith('.csv'):
                continue
            

            target_name = file.replace('.csv', '')

            #if target_name not in targets_filt:
            #    continue

            # Load the CSV file
            df = pd.read_csv(os.path.join(base_path, file))
            
            # Ensure the 'tmscore' column is in the dataframe
            if metric not in df.columns:
                print(f"Warning: {metric} column not found in file {file}. Skipping.")
                continue

            # Extract group names
            df['group'] = df['Model'].apply(extract_group_name)
            
            # Filter for Model 1 submissions
            df = df[df['Model'].str.contains('_1')]
            # print(df) 
            # Add target name to the dataframe
            
            target_name = file.replace('.csv', '')
            df['target'] = target_name
             
            # Calculate initial z-scores
            df[f"Z_{metric}"] = calculate_z_scores(df, metric)
            
            # Remove outliers (Z-scores below tolerance threshold)
            df = df[df[f"Z_{metric}"] >= -2.0]
            
            # Recalculate z-scores on the reduced dataset
            df[f"Z_{metric}_recalculated"] = calculate_z_scores(df, metric)
            
            # Apply penalty threshold
            df[f"Z_{metric}_final"] = df[f"Z_{metric}_recalculated"]#.apply(
            #    lambda x: max(x, 0)
            #)
            
            # Retain necessary columns
            metric_results.append(df[['target', 'group', f"Z_{metric}_final"]])

        
        # Combine all results for the metric
        if metric_results:
            combined_df = pd.concat(metric_results, ignore_index=True)
            
            # Pivot table: rows as targets, columns as groups
            pivot_df = combined_df.pivot(index='target', columns='group', values=f"Z_{metric}_final")
            
            pivot_df.fillna(0, inplace=True)

            # Save the pivoted table to a CSV file
            output_file = f"per_target_z_scores_{metric}_official.csv"
            pivot_df.to_csv(output_file)
            print(f"Per-target Z-score table for metric '{metric}' saved to '{output_file}'.")

# Base path to the metrics folders
base_path = "converted_csvs/"

# Metrics folders
metrics = ["GDT_TS"]

# Generate per-target Z-score tables
generate_per_target_z_score_tables(base_path, metrics)
