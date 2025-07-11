#!/bin/bash
#SBATCH --job-name=plasmids_split
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4 
#SBATCH --mem=8G  
#SBATCH --partition=allgroups
#SBATCH --time=10-00  


cd /nfsd/bcb/bcbg/zareshahraki/scripts/

/home/zareshahra/miniconda3/envs/shabnam/bin/python3.12 plasmids_split_set.py


echo "Split done!"
