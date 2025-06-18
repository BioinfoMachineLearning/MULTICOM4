#!/bin/bash

# Script to automate the monomer prediction pipeline
# Author: Your Name
# Usage: ./run_monomer_pipeline.sh <option_file> <fasta_path> <output_dir>
# Example: ./run_monomer_pipeline.sh bin/db_option_default_local /path/to/T1200.fasta /path/to/output_dir

# Check for input arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <option_file> <fasta_path> <output_dir>"
    echo "Example: $0 bin/db_option_default_local /bmlfast/casp16_ts_qs_qa/TS_run/fasta/T1200.fasta /bmlfast/casp16_ts_qs_qa/TS_run/valid/T1200"
    exit 1
fi

# Set variables from input arguments
OPTION_FILE=$1
FASTA_PATH=$2
OUTPUT_DIR=$3

# Extract the target ID from the fasta file name
TARGET_ID=$(basename "$FASTA_PATH" .fasta)

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Function to check the exit status of the previous command
check_status() {
    if [ $? -ne 0 ]; then
        echo "Error occurred in step: $1. Exiting."
        exit 1
    fi
}

# Step 1: Generate alignment
echo "Step 1: Generating alignment for target ${TARGET_ID}..."
python bin/monomer/monomer_alignment.py --option_file "$OPTION_FILE" --fasta_path "$FASTA_PATH" --output_dir "$OUTPUT_DIR"
check_status "monomer_alignment.py"

# Step 2: Generate AlphaFold2 predictors
echo "Step 2: Running AlphaFold2 predictors for target ${TARGET_ID}..."
python bin/monomer/monomer_af.py --option_file "$OPTION_FILE" --fasta_path "$FASTA_PATH" --output_dir "$OUTPUT_DIR"
check_status "monomer_af.py"

# Step 3: Run non-AlphaFold2 predictors
# echo "Step 3: Running non-AlphaFold2 predictors for target ${TARGET_ID}..."
# python bin/monomer/monomer_non_af.py --option_file "$OPTION_FILE" --fasta_path "$FASTA_PATH" --output_dir "$OUTPUT_DIR"
# check_status "monomer_non_af.py"

# Step 4: Run post-processing scripts
echo "Step 4: Running post-processing after AlphaFold2 for target ${TARGET_ID}..."
python bin/monomer/monomer_post_af.py --option_file "$OPTION_FILE" --fasta_path "$FASTA_PATH" --output_dir "$OUTPUT_DIR"
check_status "monomer_post_af.py"

# Step 5: Rank models based on finished predictors
echo "Step 5: Ranking models for target ${TARGET_ID}..."
python bin/monomer/monomer_ranking.py --option_file "$OPTION_FILE" --fasta_path "$FASTA_PATH" --output_dir "$OUTPUT_DIR"
check_status "monomer_ranking.py"

# Final message
echo "Pipeline completed successfully for target ${TARGET_ID}."
echo "Outputs are available in: ${OUTPUT_DIR}"

