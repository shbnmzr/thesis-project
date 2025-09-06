#!/usr/bin/env python3
"""
Split raw genome files into train/val/test sets.

Usage:
    python split_sets.py \
        --category plasmid \
        --base_dir ./dataset \
        --train_ratio 0.85 \
        --val_ratio 0.05 \
        --test_ratio 0.10
"""

import argparse
import random
import shutil
from pathlib import Path


def split_files(category: str, base_dir: Path, train_ratio: float, val_ratio: float, test_ratio: float, seed: int = 42):
    """
    Split genome files into train/val/test folders.

    Args:
        category (str): Dataset category (e.g., plasmid, archaea).
        base_dir (Path): Base dataset directory.
        train_ratio (float): Proportion of training data.
        val_ratio (float): Proportion of validation data.
        test_ratio (float): Proportion of test data.
        seed (int): Random seed for reproducibility.
    """
    raw_dir = base_dir / "train" / category / "raw"
    if not raw_dir.exists():
        raise FileNotFoundError(f"Raw directory not found: {raw_dir}")

    random.seed(seed)
    all_files = sorted(raw_dir.glob("*.fna"))
    all_files = list(all_files)
    random.shuffle(all_files)

    total = len(all_files)
    cut1 = int(total * train_ratio)
    cut2 = cut1 + int(total * val_ratio)

    split_dirs = {
        "train": all_files[:cut1],
        "val": all_files[cut1:cut2],
        "test": all_files[cut2:]
    }

    for split, files in split_dirs.items():
        target_dir = base_dir / split / category
        target_dir.mkdir(parents=True, exist_ok=True)
        for file in files:
            shutil.copy(file, target_dir / file.name)

    print(f"[DONE] Split {total} '{category}' genomes into train/val/test.")
    print(f"  Train: {len(split_dirs['train'])} | Val: {len(split_dirs['val'])} | Test: {len(split_dirs['test'])}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split genome files into train/val/test sets.")
    parser.add_argument("--category", required=True, help="Dataset category (e.g., plasmid, archaea).")
    parser.add_argument("--base_dir", default="./4CAC_dataset", help="Base dataset directory.")
    parser.add_argument("--train_ratio", type=float, default=0.85, help="Proportion for training set.")
    parser.add_argument("--val_ratio", type=float, default=0.05, help="Proportion for validation set.")
    parser.add_argument("--test_ratio", type=float, default=0.10, help="Proportion for test set.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility.")
    args = parser.parse_args()

    # Check that ratios sum to ~1
    total_ratio = args.train_ratio + args.val_ratio + args.test_ratio
    if abs(total_ratio - 1.0) > 1e-6:
        raise ValueError(f"Ratios must sum to 1.0 (got {total_ratio}).")

    split_files(
        args.category,
        Path(args.base_dir),
        args.train_ratio,
        args.val_ratio,
        args.test_ratio,
        args.seed
    )
