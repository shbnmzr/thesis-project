#!/bin/bash
#SBATCH --job-name=protozoa_batch
#SBATCH --output=logs/output_protozoa_%A_%a.log
#SBATCH --error=logs/error_protozoa_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=allgroups
#SBATCH --array=0-15  # Adjust based on total protozoa genome count (e.g., 1608 / 100 = 17 → 0-16)
#SBATCH --time=1-00:00:00

cd /nfsd/bcb/bcbg/zareshahraki/scripts/

START=$((SLURM_ARRAY_TASK_ID * 100))
END=$((START + 100))

/home/zareshahra/miniconda3/envs/shabnam/bin/python3.12 protozoa_process_batch.py --start $START --end $END
