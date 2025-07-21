#!/usr/bin/env python3
"""
Metagenomic Assembly Pipeline for Long Reads
--------------------------------------------
This script assembles Nanopore long reads using metaFlye.

Features:
- Logs all steps to console and log file.
- Checks for existing outputs to avoid redundant runs.
- Captures stdout and stderr for reproducibility.
- Explicit error handling and validation of inputs.

Example:
    python assemble_long_reads_metaflye.py \
        --reads ./simulated_long/generic_nanopore.fastq \
        --output_dir ./assembly_generic_flye \
        --threads 32
"""

import argparse
import subprocess
import logging
from pathlib import Path
import sys
from scripts.utils.logger import setup_logger


def check_existing_outputs(output_dir: Path):
    """Check if metaFlye assembly has already been completed."""
    assembly = output_dir / "assembly.fasta"
    graph = output_dir / "assembly_graph.gfa"
    if assembly.exists() and graph.exists():
        logging.info(f"Assembly already completed: {assembly} and {graph} exist.")
        return True
    return False


def run_metaflye(reads: Path, output_dir: Path, threads: int):
    """Run metaFlye assembler."""
    cmd = [
        "flye",
        "--nano-raw", str(reads),
        "--meta",
        "--out-dir", str(output_dir),
        "--threads", str(threads)
    ]

    logging.info(f"Running metaFlye with command:\n{' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logging.info(result.stdout)
        logging.info("metaFlye completed successfully.")
        logging.info(f"Output files: {output_dir / 'assembly.fasta'}, {output_dir / 'assembly_graph.gfa'}")
    except subprocess.CalledProcessError as e:
        logging.error("metaFlye failed!")
        logging.error(f"Return Code: {e.returncode}")
        logging.error(f"Standard Output:\n{e.stdout}")
        logging.error(f"Standard Error:\n{e.stderr}")
        sys.exit(1)


def validate_inputs(reads: Path):
    """Validate that input FASTQ file exists and is not empty."""
    if not reads.exists() or reads.stat().st_size == 0:
        logging.error(f"Input FASTQ file missing or empty: {reads}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Assemble metagenomic long reads using metaFlye.")
    parser.add_argument("--reads", required=True, help="Path to Nanopore long-read FASTQ file.")
    parser.add_argument("--output_dir", required=True, help="Directory for metaFlye output.")
    parser.add_argument("--threads", type=int, default=32, help="Number of CPU threads (default: 32).")
    parser.add_argument("--log_file", default="assemble_long_reads.log", help="Path to the log file.")
    args = parser.parse_args()

    setup_logger(Path(args.log_file))

    reads = Path(args.reads).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info(f"Starting metaFlye assembly:")
    logging.info(f"Reads={reads}\nOutput={output_dir}")

    validate_inputs(reads)

    if not check_existing_outputs(output_dir):
        run_metaflye(reads, output_dir, args.threads)

    logging.info("Long-read assembly pipeline finished successfully.")


if __name__ == "__main__":
    main()
