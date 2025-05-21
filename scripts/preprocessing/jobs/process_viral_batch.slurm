#!/bin/bash
#SBATCH --job-name=viral_batch
#SBATCH --output=logs/output_viral_%A_%a.log
#SBATCH --error=logs/error_viral_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --partition=allgroups
#SBATCH --array=0-213  # Adjust depending on number of 1000-sized batches
#SBATCH --time=10-00:00:00

cd /nfsd/bcb/bcbg/zareshahraki/scripts/
START=$((SLURM_ARRAY_TASK_ID * 1000))
END=$((START + 1000))

/home/zareshahra/miniconda3/envs/shabnam/bin/python3.12 viral_process_batch.py --start $START --end $END
