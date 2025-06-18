#!/bin/bash
#SBATCH --job-name=feature_archaea_test
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=allgroups
#SBATCH --time=10-00:00:00

cd /path/to/scripts/

/path/to/venv/bin/python3.12 feature_extraction.py \
  --dataset_root /path/to/4CAC_dataset \
  --log_root /path/to/logs \
  --output_dir /path/to/processed_data_jellyfish \
  --split val \
  --category archaea
echo "Feature extraction for split done!"
