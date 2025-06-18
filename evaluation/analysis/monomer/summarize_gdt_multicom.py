import os
import pandas as pd

# Define path to CSV files and group ID
base_path = "./converted_csvs"  # Adjust as needed
group_id = "TS051"
valid_suffixes = [f"{group_id}_{i}" for i in range(1, 6)]

# Store results
results = []

# Process all relevant CSV files
for file in os.listdir(base_path):
    if not file.endswith(".csv") or not file.startswith("T"):
        continue

    target_name = file.replace(".csv", "")
    file_path = os.path.join(base_path, file)
    try:
        df = pd.read_csv(file_path)

        if "Model" in df.columns and "GDT_TS" in df.columns:
            filtered = df[df["Model"].str.contains(group_id) & df["Model"].str.contains("_")]
            filtered = filtered[filtered["Model"].apply(lambda x: any(f in x for f in valid_suffixes))]

            top1_row = filtered[filtered["Model"].str.contains(f"{group_id}_1")]
            top1_tmscore = top1_row["GDT_TS"].values[0] if not top1_row.empty else None
            max_tmscore = filtered["GDT_TS"].max() if not filtered.empty else None

            print(f"{target_name}\t{top1_tmscore}\t{max_tmscore}")
    except Exception as e:
        print(f"Error processing {file}: {e}")


print(results)
