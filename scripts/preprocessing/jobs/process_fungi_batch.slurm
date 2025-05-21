#!/bin/bash
#SBATCH --job-name=fungi_batch
#SBATCH --output=../4CAC_dataset/logs/output_fungi_%A_%a.log
#SBATCH --error=../4CAC_dataset/logs/error_fungi_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16                
#SBATCH --mem=32G                          # Higher memory for large genomes and CD-HIT
#SBATCH --partition=allgroups
#SBATCH --array=0-10                       # Adjust depending on total fungi files / 1000
#SBATCH --time=3-00:00:00

# Activate environment if needed
# source ~/.bashrc
# conda activate shabnam

cd /nfsd/bcb/bcbg/zareshahraki/scripts/

START=$((SLURM_ARRAY_TASK_ID * 1000))
END=$((START + 1000))

/home/zareshahra/miniconda3/envs/shabnam/bin/python3.12 fungi_process_batch.py --start $START --end $END
