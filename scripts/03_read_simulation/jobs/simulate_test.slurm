#!/bin/bash
#SBATCH --job-name=simulate_test
#SBATCH --output=logs/slurm/simulate_%x_%j.out
#SBATCH --error=logs/slurm/simulate_%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=10-00:00:00
#SBATCH --partition=allgroups
#SBATCH --tmp=8G

# === Load Environment ===
set -euo pipefail
echo "[START] Job ID: $SLURM_JOB_ID | Host: $(hostname)"
echo "[INFO] Starting read simulation..."

# === Parameters ===
SPLIT="test"
SCRIPT_DIR="path/to/scripts"
ENV_PATH="path/to/bin/python3.9"
SCRIPT_NAME="simulate_reads.py"
GENOMES_FASTA="path/to/combined/${SPLIT}_combined.fasta"
OUTPUT_DIR="path/to/reads/${SPLIT}"
LOG_FILE="logs/simulate_${SPLIT}.log"
OUTPUT_PREFIX="${SPLIT}_sim"
N_READS=10000000

# === Directory Preparation ===
mkdir -p logs/slurm "$OUTPUT_DIR"

# === Launch Python Read Simulation Script ===
cd "$SCRIPT_DIR"
"$ENV_PATH" "$SCRIPT_NAME" \
    --genomes_fasta "$GENOMES_FASTA" \
    --output_dir "$OUTPUT_DIR" \
    --output_prefix "$OUTPUT_PREFIX" \
    --n_reads "$N_READS" \
    --cpus 8 \
    --log_file "$LOG_FILE"

echo "[DONE] Completed read simulation for $SPLIT at $(date)"
