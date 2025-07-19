#!/usr/bin/env python3
"""
Combine all cleaned FASTA files from multiple categories into a single FASTA file.

Expected directory structure:
    dataset_dir/
        archaea/cleaned/*.fna
        bacteria/cleaned/*.fna
        ...

Usage:
    python combine_fasta.py \
        --dataset_dir ./4CAC_dataset/train \
        --output_file ./4CAC_dataset/combined/all_cleaned.fasta \
        --log_file ./logs/combine_fasta.log
"""

import argparse
import logging
from pathlib import Path
from scripts.utils.logger import setup_logger


def combine_fasta(split_dir: Path, output_file: Path, categories=None):
    """
    Combine all .fna files from cleaned directories into a single FASTA.

    Args:
        split_dir (Path): Root dataset directory (e.g., ./4CAC_dataset/train).
        output_file (Path): Output combined FASTA path.
        categories (list): List of categories to process. Defaults to standard categories.
    """
    if categories is None:
        categories = ["archaea", "bacteria", "fungi", "protozoa", "plasmid", "viral"]

    fasta_files = []
    for category in categories:
        cleaned_dir = split_dir / category / "cleaned"
        if cleaned_dir.exists():
            found_files = list(cleaned_dir.glob("*.fna"))
            logging.info(f"Found {len(found_files)} files in {cleaned_dir}")
            fasta_files.extend(found_files)
        else:
            logging.warning(f"Directory {cleaned_dir} does not exist. Skipping.")

    if not fasta_files:
        logging.error(f"No .fna files found under {split_dir}. Exiting.")
        raise FileNotFoundError(f"No .fna files found under {split_dir}")

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as outfile:
        for fasta in sorted(fasta_files):
            logging.info(f"Adding {fasta.name}")
            with open(fasta, "r") as infile:
                outfile.write(infile.read())

    logging.info(f"✅ Combined {len(fasta_files)} FASTA files into {output_file.resolve()}")


def main():
    parser = argparse.ArgumentParser(description="Combine cleaned FASTA files into a single FASTA.")
    parser.add_argument("--dataset_dir", required=True, help="Root dataset directory (e.g., ./4CAC_dataset/train)")
    parser.add_argument("--output_file", required=True, help="Path to save the combined FASTA.")
    parser.add_argument("--log_file", default="./logs/combine_fasta.log", help="Path to the log file.")
    parser.add_argument("--categories", nargs="+", default=None,
                        help="Optional: List of categories to include (default: standard categories).")
    args = parser.parse_args()

    setup_logger(Path(args.log_file))

    split_dir = Path(args.dataset_dir).resolve()
    output_file = Path(args.output_file).resolve()

    logging.info(f"Starting FASTA combination from {split_dir} to {output_file}")
    combine_fasta(split_dir, output_file, args.categories)
    logging.info("Done.")


if __name__ == "__main__":
    main()
