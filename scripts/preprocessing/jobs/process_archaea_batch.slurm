#!/bin/bash
#SBATCH --job-name=archaea_batch
#SBATCH --output=../4CAC_dataset/logs/output_archaea_%A_%a.log
#SBATCH --error=../4CAC_dataset/logs/error_archaea_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=allgroups
#SBATCH --array=0-12               # Adjust based on total archaea files / 1000
#SBATCH --time=2-00:00:00

# Activate environment if needed
# source ~/.bashrc
# conda activate shabnam

cd /nfsd/bcb/bcbg/zareshahraki/scripts/

START=$((SLURM_ARRAY_TASK_ID * 1000))
END=$((START + 1000))

/home/zareshahra/miniconda3/envs/shabnam/bin/python3.12 archaea_process_batch.py --start $START --end $END
