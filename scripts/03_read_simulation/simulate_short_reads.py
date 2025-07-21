#!/usr/bin/env python3
"""
Simulate Illumina short reads

Features:
1. Uses InSilicoSeq to generate realistic metagenomic Illumina reads.
2. Supports multiple scenarios (generic + filtered).
3. Outputs standard FASTQ files for downstream assembly, testing, or GNN training.
4. Saves metadata for reproducibility.

Example:
    python simulate_short_reads.py \
        --fasta ./4CAC_dataset/test_prep/all_references.fasta \
        --abundance_generic ./4CAC_dataset/test_prep/abundance_generic.txt \
        --abundance_filtered ./4CAC_dataset/test_prep/abundance_filtered.txt \
        --output_dir ./4CAC_dataset/simulated_short \
        --model miseq --n_reads 10000000 --seed 42
"""

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path
from scripts.utils.logger import setup_logger

# -----------------------------
# ISS WRAPPER
# -----------------------------
def run_iss(fasta: Path, abundance: Path, output_prefix: Path,
            n_reads: int, cpus: int, model: str):
    """
    Run InSilicoSeq to simulate paired-end Illumina reads.

    Args:
        fasta (Path): Combined reference FASTA file.
        abundance (Path): Abundance file for this scenario.
        output_prefix (Path): Prefix for output FASTQs.
        n_reads (int): Total number of reads to simulate.
        cpus (int): CPU threads to use.
        model (str): Illumina model (miseq, hiseq, nova, etc.).
    """
    cmd = [
        "iss", "generate",
        "--genomes", str(fasta),
        "--abundance_file", str(abundance),
        "--model", model,
        "--output", str(output_prefix),
        "--n_reads", str(n_reads),
        "--cpus", str(cpus)
    ]
    logging.info(f"Running InSilicoSeq: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
        logging.info(f"Simulation complete: {output_prefix}_R1.fastq, {output_prefix}_R2.fastq")
    except subprocess.CalledProcessError as e:
        logging.error(f"InSilicoSeq failed: {e}")
        sys.exit(1)
    except FileNotFoundError:
        logging.error("InSilicoSeq (iss) is not installed or not in PATH.")
        sys.exit(1)

# -----------------------------
# MAIN
# -----------------------------
def main():
    parser = argparse.ArgumentParser(description="Simulate Illumina short reads for metagenomic test sets.")
    parser.add_argument("--fasta", required=True, help="Combined reference FASTA file (with class prefixes).")
    parser.add_argument("--abundance_generic", required=True, help="Abundance file for generic scenario.")
    parser.add_argument("--abundance_filtered", required=True, help="Abundance file for filtered scenario.")
    parser.add_argument("--output_dir", required=True, help="Directory to save simulated reads.")
    parser.add_argument("--model", default="miseq", help="Illumina model (miseq, hiseq, nova, etc.).")
    parser.add_argument("--n_reads", type=int, default=10000000, help="Total reads per scenario.")
    parser.add_argument("--cpus", type=int, default=8, help="Number of CPU threads to use.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility.")
    parser.add_argument("--log_file", default="./logs/simulate_short_reads.log",
                        help="Log file (default: ./logs/simulate_short_reads.log)")
    args = parser.parse_args()

    # Setup logging
    setup_logger(Path(args.log_file))

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Run simulations for generic and filtered
    run_iss(Path(args.fasta), Path(args.abundance_generic),
            output_dir / "generic_miseq", args.n_reads, args.cpus, args.model)

    run_iss(Path(args.fasta), Path(args.abundance_filtered),
            output_dir / "filtered_miseq", args.n_reads, args.cpus, args.model)

    # Save metadata
    metadata = {
        "fasta": str(args.fasta),
        "abundance_generic": str(args.abundance_generic),
        "abundance_filtered": str(args.abundance_filtered),
        "model": args.model,
        "n_reads": args.n_reads,
        "cpus": args.cpus,
        "seed": args.seed
    }
    with open(output_dir / "simulation_metadata.json", "w") as f:
        json.dump(metadata, f, indent=4)

    logging.info(f"Metadata saved: {output_dir / 'simulation_metadata.json'}")
    logging.info("Short-read simulation complete for both scenarios.")


if __name__ == "__main__":
    main()
