import pandas as pd
import os

# Function to generate the final table
def generate_final_table_with_threshold(base_path, metrics, threshold=-2.0):
    final_results = []
    
    for metric in metrics:
        # File path for the per-target Z-score table
        input_file = f"per_target_z_scores_{metric}_official.csv"
        input_path = os.path.join(base_path, input_file)
        
        if not os.path.exists(input_path):
            print(f"Warning: File '{input_file}' does not exist. Skipping.")
            continue
        
        # Load the per-target Z-score table
        df = pd.read_csv(input_path, index_col=0)
        
        # Sum Z-scores for each group across all targets
        group_sum = df.sum(axis=0).reset_index()
        group_sum.columns = ["group", f"Sum_Z_{metric}"]
        
        final_results.append(group_sum)
    
    # Merge results for all metrics
    final_table = pd.DataFrame()
    for group_sum in final_results:
        if final_table.empty:
            final_table = group_sum
        else:
            final_table = pd.merge(final_table, group_sum, on="group", how="outer")
    
    # Fill NaN values with 0 for groups missing data in some metrics
    final_table.fillna(0, inplace=True)
    
    # Add a column for the total sum of Z-scores across all metrics
    z_score_columns = [col for col in final_table.columns if col.startswith("Sum_Z_")]
    final_table["Total_Sum_Z_Scores"] = final_table[z_score_columns].sum(axis=1)
    
    # Save the final table to a CSV file
    output_file = os.path.join(base_path, "final_z_scores_table_with_threshold_official.csv")
    final_table.to_csv(output_file, index=False)
    print(f"Final Z-score table with threshold saved to '{output_file}'.")

# Base path where per-target Z-score tables are stored
base_path = "."

# Metrics
metrics = ["GDT_TS"]

# Generate the final table with threshold
generate_final_table_with_threshold(base_path, metrics, threshold=0)
