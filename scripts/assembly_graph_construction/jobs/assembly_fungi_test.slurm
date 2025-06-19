#!/bin/bash
#SBATCH --job-name=spades_fungi_test
#SBATCH --output=spades_%x_%j.out
#SBATCH --error=spades_%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=10-00:00:00
#SBATCH --partition=allgroups

# === PARAMETERS (edit these) ===
SPLIT="test"
CATEGORY="fungi"

# Run the updated Python script
cd /nfsd/bcb/bcbg/zareshahraki/scripts/

/home/zareshahra/miniconda3/envs/shabnam/bin/python3.12 assembly.py \
    --split "$SPLIT" \
    --category "$CATEGORY"