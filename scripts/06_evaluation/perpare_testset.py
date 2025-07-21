#!/usr/bin/env python3
"""
Prepare metagenomic test set for benchmarking classification tools.

Features:
1. Combine genomes from class-organized directories into a single FASTA (with class prefixes).
2. Generate abundance profiles for different scenarios (4CAC-style by default).
3. Output metadata for reproducibility.

This script is generic: you can specify custom class ratios and distributions.

Example (4CAC-style defaults):
    python prepare_testset.py \
        --base_dir ./datasets/genomes \
        --output_dir ./datasets/test_prep \
        --preset 4cac \
        --seed 42
"""

import argparse
import logging
import json
import numpy as np
from math import log10
from pathlib import Path
from typing import Dict, List, Tuple
import random
from scripts.utils.logger import setup_logger

# -----------------------------
# DEFAULT PRESETS
# -----------------------------
PRESETS = {
    "4cac": {
        "generic_ratio": {"prokaryote": 56, "microeuk": 24, "virus": 10, "plasmid": 10},
        "filtered_ratio": {"prokaryote": 14, "microeuk": 6, "virus": 40, "plasmid": 40},
        "distributions": {"default": "lognormal", "plasmid": "geometric"}
    },
    "uniform": {
        "generic_ratio": {"prokaryote": 25, "microeuk": 25, "virus": 25, "plasmid": 25},
        "filtered_ratio": {"prokaryote": 25, "microeuk": 25, "virus": 25, "plasmid": 25},
        "distributions": {"default": "uniform", "plasmid": "uniform"}
    }
}

# -----------------------------
# CLASS MAP
# -----------------------------
CLASS_MAP = {
    "bacteria": "prokaryote",
    "archaea": "prokaryote",
    "viral": "virus",
    "fungi": "microeuk",
    "protozoa": "microeuk",
    "plasmid": "plasmid"
}

# -----------------------------
# DISTRIBUTION FUNCTIONS
# -----------------------------
def generate_distribution(dist_type: str, n: int) -> np.ndarray:
    """Generate a normalized abundance profile."""
    if dist_type == "lognormal":
        vals = np.random.lognormal(mean=1.0, sigma=2.0, size=n)
    elif dist_type == "uniform":
        vals = np.ones(n)
    else:
        raise ValueError(f"Unsupported distribution type: {dist_type}")
    return vals / vals.sum()

def geometric_copy_number(length_bp: int) -> int:
    """Generate plasmid copy number (4CAC logic)."""
    p = min(1, log10(length_bp) / 7) if length_bp > 0 else 1
    return np.random.geometric(p)

# -----------------------------
# CORE LOGIC
# -----------------------------
def combine_genomes(base_dir: Path, output_fasta: Path) -> Dict[str, Tuple[str, int]]:
    """
    Combine genome FASTAs into one with class prefixes.

    Returns:
        dict: {seq_id: (class_name, sequence_length)}
    """
    seq_lengths = {}
    with open(output_fasta, "w") as out:
        for folder, cls in CLASS_MAP.items():
            genome_dir = base_dir / folder
            if not genome_dir.exists():
                logging.warning(f"Skipping missing folder: {genome_dir}")
                continue
            for fasta in genome_dir.glob("*.fna"):
                with open(fasta) as f:
                    seq_id, seq_data = None, []
                    for line in f:
                        if line.startswith(">"):
                            if seq_id:
                                seq_lengths[seq_id] = (cls, len("".join(seq_data)))
                                out.write(f">{seq_id}\n{''.join(seq_data)}\n")
                            raw_id = line.strip()[1:].split()[0]
                            seq_id = f"{cls}|{raw_id}"
                            seq_data = []
                        else:
                            seq_data.append(line.strip())
                    if seq_id:
                        seq_lengths[seq_id] = (cls, len("".join(seq_data)))
                        out.write(f">{seq_id}\n{''.join(seq_data)}\n")
    logging.info(f"Combined FASTA created: {output_fasta}")
    return seq_lengths

def generate_abundances(seq_lengths: Dict[str, Tuple[str, int]],
                        ratios: Dict[str, float],
                        distributions: Dict[str, str]) -> Dict[str, float]:
    """
    Generate normalized abundances for all sequences based on class ratios.
    """
    abundances = {}
    # Non-plasmid classes
    for cls, target_ratio in ratios.items():
        if cls == "plasmid":
            continue
        seqs = [s for s, (c, _) in seq_lengths.items() if c == cls]
        if not seqs:
            continue
        dist = generate_distribution(distributions.get("default", "lognormal"), len(seqs))
        for s, a in zip(seqs, dist):
            abundances[s] = a * target_ratio

    # Plasmids
    plasmids = [s for s, (c, _) in seq_lengths.items() if c == "plasmid"]
    if plasmids:
        if distributions.get("plasmid") == "geometric":
            # Approximate host abundance per plasmid (prokaryote share / plasmid count)
            host_share = ratios.get("prokaryote", 1) / max(len(plasmids), 1)
            for s in plasmids:
                _, length = seq_lengths[s]
                abundances[s] = host_share * geometric_copy_number(length)
        else:  # uniform plasmids
            for s in plasmids:
                abundances[s] = ratios.get("plasmid", 1) / len(plasmids)

    # Normalize to sum 1
    total = sum(abundances.values())
    for s in abundances:
        abundances[s] /= total
    return abundances

def save_abundance_file(abundance: Dict[str, float], out_file: Path):
    with open(out_file, "w") as f:
        for s, a in abundance.items():
            f.write(f"{s}\t{a:.8f}\n")
    logging.info(f"Abundance file saved: {out_file}")

# -----------------------------
# MAIN
# -----------------------------
def main():
    parser = argparse.ArgumentParser(description="Prepare metagenomic test set (generic & filtered scenarios)")
    parser.add_argument("--base_dir", required=True, help="Base directory with class-organized genome folders")
    parser.add_argument("--output_dir", required=True, help="Output directory for FASTA and abundance files")
    parser.add_argument("--preset", choices=PRESETS.keys(), default="4cac", help="Preset scenario (default: 4cac)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    parser.add_argument("--log_file", default="./logs/prepare_testset.log",
                        help="Path to log file (default: ./logs/prepare_testset.log)")

    args = parser.parse_args()

    setup_logger(Path(args.log_file))

    # Set seeds
    random.seed(args.seed)
    np.random.seed(args.seed)

    base_dir = Path(args.base_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    combined_fasta = output_dir / "all_references.fasta"
    seq_lengths = combine_genomes(base_dir, combined_fasta)

    preset = PRESETS[args.preset]
    # Generic
    generic_ab = generate_abundances(seq_lengths, preset["generic_ratio"], preset["distributions"])
    save_abundance_file(generic_ab, output_dir / "abundance_generic.txt")
    # Filtered
    filtered_ab = generate_abundances(seq_lengths, preset["filtered_ratio"], preset["distributions"])
    save_abundance_file(filtered_ab, output_dir / "abundance_filtered.txt")

    # Save metadata
    meta = {
        "preset": args.preset,
        "seed": args.seed,
        "genome_counts": {cls: sum(1 for c in seq_lengths.values() if c[0] == cls) for cls in CLASS_MAP.values()}
    }
    with open(output_dir / "metadata.json", "w") as f:
        json.dump(meta, f, indent=4)
    logging.info(f"Metadata saved: {output_dir / 'metadata.json'}")
    logging.info("Test set preparation complete.")


if __name__ == "__main__":
    main()
