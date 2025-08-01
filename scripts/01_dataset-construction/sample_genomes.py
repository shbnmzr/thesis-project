#!/usr/bin/env python3
"""
Randomly samples N (a given number) genomes from a directory.

Usage:
    python sample_genomes.py --source_dir /path/to/train/bacteria --output_dir ./my_dataset --n 1000
"""

import argparse
import random
from pathlib import Path
import shutil

def sample_files(source_dir: Path, output_dir: Path, n: int, seed: int):
    """
    Samples n genomes from a directory
    source_dir: source directory
    output_dir: output directory
    n: number of genomes to sample
    seed: random seed
    """
    fasta_files = list(source_dir.glob("*.fna")) + list(source_dir.glob("*.fa")) + list(source_dir.glob("*.fasta"))

    if not fasta_files:
        print(f"No FASTA files found in {source_dir}")
        return

    if n > len(fasta_files):
        raise ValueError(f"Requested {n} files, but only found {len(fasta_files)} in {source_dir}")

    random.seed(seed)
    sampled = random.sample(fasta_files, n)

    output_dir.mkdir(parents=True, exist_ok=True)
    for f in sampled:
        shutil.copy(f, output_dir / f.name)

    with open(output_dir / "sampled_filenames.txt", "w") as log:
        for f in sampled:
            log.write(f.name + "\n")

    print(f"Sampled {n} genomes from {source_dir.name}. Output: {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="Sample genome files from local folders.")
    parser.add_argument("--source_dir", required=True, help="Directory with downloaded genomes (e.g., train/bacteria)")
    parser.add_argument("--output_dir", required=True, help="Directory to save sampled genomes")
    parser.add_argument("--n", type=int, required=True, help="Number of genomes to sample")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    sample_files(Path(args.source_dir), Path(args.output_dir), args.n, args.seed)

if __name__ == "__main__":
    main()
