import os
import random
from pathlib import Path
import shutil

random.seed(42)

RAW_DIR = Path("../4CAC_dataset/train/plasmid/raw")
BASE_DIR = Path("../4CAC_dataset")
SPLITS = {
    "train": 0.85,
    "val": 0.05,
    "test": 0.10
}

all_files = sorted(list(RAW_DIR.glob("*.fna")))
random.shuffle(all_files)

total = len(all_files)
cut1 = int(total * SPLITS["train"])
cut2 = cut1 + int(total * SPLITS["val"])

split_dirs = {
    "train": all_files[:cut1],
    "val": all_files[cut1:cut2],
    "test": all_files[cut2:]
}

for split, files in split_dirs.items():
    target_dir = BASE_DIR / split / "plasmid"
    target_dir.mkdir(parents=True, exist_ok=True)
    for file in files:
        shutil.copy(file, target_dir / file.name)

print(f"[DONE] Split {total} plasmid genomes into train/val/test.")