import pandas as pd
import re
from io import StringIO
import os

outdir = 'converted_csvs'
os.makedirs(outdir, exist_ok=True)

for target_result_file in os.listdir('downloaded_files'):
    print("processing " + target_result_file)
    # Simulated user input string (shortened here for practicality; assume full string in actual code)
    lines = open('downloaded_files/' + target_result_file).readlines()
    if len(lines) <= 5:
        continue
    # Remove the comment sign and extra whitespace from the header
    header = lines[0].split()
    data_lines = lines[1:]
    csv_data = ''.join(data_lines)

    # Convert to DataFrame
    df = pd.read_csv(StringIO(csv_data), delim_whitespace=True, header=None)
    df.columns = header
    
    df.to_csv(outdir + '/' + target_result_file.replace('.txt', '.csv'))




