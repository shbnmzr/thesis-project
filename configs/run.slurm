#!/bin/bash
#SBATCH --job-name thesis_zareshahraki
#SBATCH --output output_%j.log
#SBATCH --error errors_%j.log
#SBATCH --ntasks 1
#SBATCH --mem 16G
#SBATCH --partition allgroups
#SBATCH --gres=gpu:rtx
#SBATCH --time 5-00

cd /home/zareshara/thesis
srun /home/zareshahra/miniconda3/envs/shabnam/bin/python3.12 app.py
