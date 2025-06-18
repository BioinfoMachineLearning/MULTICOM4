#!/bin/bash

# Script to automate the multimer prediction pipeline
# Author: Your Name
# Usage: ./run_multimer_pipeline.sh <option_file> <fasta_path> <output_dir>
# Example: ./run_multimer_pipeline.sh bin/db_option_default_local /path/to/H1208.fasta /path/to/output_dir

# Check for input arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <option_file> <fasta_path> <output_dir>"
    echo "Example: $0 bin/db_option_default_local /bmlfast/casp16_ts_qs_qa/TS_run/fasta/H1208.fasta /bmlfast/casp16_ts_qs_qa/TS_run/valid/H1208"
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

# Step 1: Generate alignment for multimer target
echo "Step 1: Generating alignment for multimer target ${TARGET_ID}..."
python bin/multimer/multimer_subunit_alignment.py --option_file "$OPTION_FILE" --fasta_path "$FASTA_PATH" --output_dir "$OUTPUT_DIR"
check_status "multimer_subunit_alignment.py"

# Step 2: Run AlphaFold2 predictors
echo "Step 2: Running AlphaFold2 predictors for target ${TARGET_ID}..."
python bin/multimer/multimer_subunit_af.py --option_file "$OPTION_FILE" --fasta_path "$FASTA_PATH" --output_dir "$OUTPUT_DIR"
check_status "multimer_subunit_af.py"

# Step 3: Run post-processing for MSA from AlphaFold2
echo "Step 3: Running post-processing after AlphaFold2 for target ${TARGET_ID}..."
python bin/multimer/multimer_subunit_post_af.py --option_file "$OPTION_FILE" --fasta_path "$FASTA_PATH" --output_dir "$OUTPUT_DIR"
check_status "multimer_subunit_post_af.py"

# Step 4: Run multimer predictors requiring MSA only
echo "Step 4: Running multimer predictors requiring only MSA for target ${TARGET_ID}..."
python bin/multimer/multimer_af_local.py --option_file "$OPTION_FILE" --fasta_path "$FASTA_PATH" --output_dir "$OUTPUT_DIR"
check_status "multimer_af_local.py"

# Step 5: Generate multimer alignments and templates
echo "Step 5: Generating multimer alignments and templates for target ${TARGET_ID}..."
python bin/multimer/multimer_alignment_af.py --option_file "$OPTION_FILE" --fasta_path "$FASTA_PATH" --output_dir "$OUTPUT_DIR"
check_status "multimer_alignment_af.py"

# Step 6: Run foldseek-based predictors
echo "Step 6: Running foldseek-based predictors for target ${TARGET_ID}..."
python bin/multimer/multimer_post_af.py --option_file "$OPTION_FILE" --fasta_path "$FASTA_PATH" --output_dir "$OUTPUT_DIR"
check_status "multimer_post_af.py"

# Step 7: Rank models based on the predictors
echo "Step 7: Ranking models for target ${TARGET_ID}..."
python bin/multimer/multimer_ranking.py --option_file "$OPTION_FILE" --fasta_path "$FASTA_PATH" --output_dir "$OUTPUT_DIR"
check_status "multimer_ranking.py"

# Final message
echo "Multimer pipeline completed successfully for target ${TARGET_ID}."
echo "Outputs are available in: ${OUTPUT_DIR}"

