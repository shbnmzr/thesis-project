#!/usr/bin/env python3
"""
Simulate Illumina reads using InSilicoSeq (iss generate).

Usage:
    python simulate_reads.py \
        --genomes_fasta ./4CAC_dataset/combined/all_cleaned.fasta \
        --output_dir ./simulated_reads \
        --output_prefix simulated_miseq \
        --n_reads 10000000 \
        --cpus 8 \
        --model miseq \
        --log_file ./logs/simulate_reads.log
"""

import argparse
import subprocess
import logging
from pathlib import Path
import sys


# -----------------------------
# LOGGING
# -----------------------------
def setup_logger(log_file: Path):
    """Set up logging to both file and stdout."""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logging.info(f"Logging initialized. Log file: {log_file.resolve()}")


# -----------------------------
# MAIN FUNCTIONALITY
# -----------------------------
def simulate_reads(genomes_fasta: Path, output_prefix: Path, n_reads: int, cpus: int, model: str):
    """
    Run InSilicoSeq to simulate Illumina reads.

    Args:
        genomes_fasta (Path): Path to input FASTA.
        output_prefix (Path): Output prefix (directory/prefix for FASTQs).
        n_reads (int): Total number of reads to simulate.
        cpus (int): Number of CPU threads.
        model (str): Sequencing model (e.g., miseq, hiseq, nova).
    """
    cmd = [
        "iss", "generate",
        "--genomes", str(genomes_fasta),
        "--model", model,
        "--output", str(output_prefix),
        "--n_reads", str(n_reads),
        "--cpus", str(cpus)
    ]

    logging.info(f"Running InSilicoSeq: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
        r1_file = f"{output_prefix}_R1.fastq"
        r2_file = f"{output_prefix}_R2.fastq"
        logging.info(f"✅ Simulation complete. Reads saved at:\n  {r1_file}\n  {r2_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"❌ InSilicoSeq failed with error: {e}")
        raise
    except FileNotFoundError:
        logging.error("❌ InSilicoSeq (iss) is not installed or not in PATH.")
        raise


def main():
    parser = argparse.ArgumentParser(description="Simulate Illumina reads using InSilicoSeq.")
    parser.add_argument("--genomes_fasta", required=True, help="Path to the combined genomes FASTA file.")
    parser.add_argument("--output_dir", required=True, help="Directory to save the simulated reads.")
    parser.add_argument("--output_prefix", default="simulated", help="Prefix for the output FASTQ files.")
    parser.add_argument("--n_reads", type=int, default=10000000, help="Total number of reads to simulate.")
    parser.add_argument("--cpus", type=int, default=8, help="Number of CPU threads to use.")
    parser.add_argument("--model", default="miseq", help="Illumina model (miseq, hiseq, nova, etc.).")
    parser.add_argument("--log_file", default="./logs/simulate_reads.log", help="Path to the log file.")
    args = parser.parse_args()

    # Setup logging
    setup_logger(Path(args.log_file))

    # Validate input FASTA
    genomes_fasta = Path(args.genomes_fasta).resolve()
    if not genomes_fasta.exists():
        logging.error(f"Genome FASTA not found: {genomes_fasta}")
        sys.exit(1)

    # Prepare output
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = output_dir / args.output_prefix

    logging.info(f"Starting read simulation using model '{args.model}'")
    simulate_reads(genomes_fasta, output_prefix, args.n_reads, args.cpus, args.model)


if __name__ == "__main__":
    main()
