#!/bin/bash
#SBATCH --job-name=bacteria_process
#SBATCH --output=output_%j.log
#SBATCH --error=errors_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4 
#SBATCH --mem=16G  
#SBATCH --partition=allgroups
#SBATCH --time=25-00  

cd path/to/scripts/ || exit 1

/path/to/miniconda3/envs/venv/bin/python3.12 process.py \
  --base_dir /path/to/4CAC_dataset \
  --category bacteria \
  --min_len 1000 \
  --min_contigs 1 \
  --max_ambig 0.1 \
  --log_dir /path/to/logs

echo "All genome processing for 'bacteria' completed successfully!"
