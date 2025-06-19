#!/bin/bash
#SBATCH --job-name=spades_plasmid_val
#SBATCH --output=logs/slurm/spades_%x_%j.out
#SBATCH --error=logs/slurm/spades_%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=10-00:00:00
#SBATCH --partition=allgroups
#SBATCH --tmp=10G

# === Load Environment ===
set -euo pipefail
echo "[START] Job ID: $SLURM_JOB_ID | Host: $(hostname)"
echo "[INFO] Starting SPAdes graph construction..."

# === Parameters ===
SPLIT="val"
CATEGORY="plasmid"
SCRIPT_DIR="/path/to/scripts"
ENV_PATH="/path/to/venv/bin/python3.12"
SCRIPT_NAME="assembly_graph_construction.py"

# === Directory Preparation ===
mkdir -p logs/slurm

# === Launch Python Assembly Script ===
cd "$SCRIPT_DIR"
"$ENV_PATH" "$SCRIPT_NAME" \
    --split "$SPLIT" \
    --category "$CATEGORY"

echo "[DONE] Completed SPAdes job for $CATEGORY/$SPLIT at $(date)"
