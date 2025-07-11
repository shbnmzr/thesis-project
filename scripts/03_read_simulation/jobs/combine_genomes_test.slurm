#!/bin/bash
#SBATCH --job-name=combine_test
#SBATCH --output=logs/slurm/combine_%x_%j.out
#SBATCH --error=logs/slurm/combine_%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --partition=allgroups
#SBATCH --tmp=4G

# === Load Environment ===
set -euo pipefail
echo "[START] Job ID: $SLURM_JOB_ID | Host: $(hostname)"
echo "[INFO] Starting genome combination..."

# === Parameters ===
SPLIT="test"
SCRIPT_DIR="path/to/scripts"
ENV_PATH="path/to/bin/python3.12"
SCRIPT_NAME="combine_genomes.py"
DATASET_DIR="path/to/${SPLIT}"
OUTPUT_DIR="path/to/combined"
LOG_FILE="logs/combine_${SPLIT}.log"

# === Directory Preparation ===
mkdir -p logs/slurm "$OUTPUT_DIR"

# === Launch Python Combination Script ===
cd "$SCRIPT_DIR"
"$ENV_PATH" "$SCRIPT_NAME" \
    --dataset_dir "$DATASET_DIR" \
    --output_file "$OUTPUT_DIR/${SPLIT}_combined.fasta" \
    --log_file "$LOG_FILE"

echo "[DONE] Completed genome combination for $SPLIT at $(date)"
