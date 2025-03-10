#!/bin/bash
#SBATCH --job-name=dl_fungi
#SBATCH --output=output_%j.log
#SBATCH --error=errors_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G  
#SBATCH --partition=allgroups
#SBATCH --time=3-00  

# Navigate to thesis directory
cd /home/zareshahra/thesis || { echo "Directory /home/zareshahra/thesis does not exist"; exit 1; }

# Run each download script simultaneously in the background
bash ./scripts/download_fungi.sh

echo "All fungi genome downloads completed successfully!"
