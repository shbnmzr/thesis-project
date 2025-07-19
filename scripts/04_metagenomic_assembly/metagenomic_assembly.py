#!/usr/bin/env python3
"""
Metagenomic Assembly Pipeline
-----------------------------
This script assembles metagenomic paired-end reads using metaSPAdes.

Features:
- Logs all steps and outputs to both console and log file.
- Captures stdout and stderr for reproducibility.
- Checks for existing outputs to avoid redundant runs.
- Explicit error handling with clear messages.
"""

import argparse
import subprocess
import logging
from pathlib import Path
import sys
from scripts.utils.logger import setup_logger


def check_existing_outputs(output_dir: Path):
    """Check if assembly has already been completed."""
    scaffolds = output_dir / "scaffolds.fasta"
    graph = output_dir / "assembly_graph_with_scaffolds.gfa"
    if scaffolds.exists() and graph.exists():
        logging.info(f"Assembly already completed: {scaffolds} and {graph} exist.")
        return True
    return False


def run_metaspades(read1: Path, read2: Path, output_dir: Path, threads: int, memory: int):
    """Run metaSPAdes assembler."""
    cmd = [
        "spades.py", "--meta",
        "-1", str(read1),
        "-2", str(read2),
        "-o", str(output_dir),
        "--threads", str(threads),
        "--memory", str(memory),
    ]

    logging.info(f"Running metaSPAdes with command:\n{' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logging.info(result.stdout)
        logging.info(f"metaSPAdes completed successfully.")
        logging.info(f"Output files: {output_dir / 'scaffolds.fasta'}, {output_dir / 'assembly_graph_with_scaffolds.gfa'}")
    except subprocess.CalledProcessError as e:
        logging.error("metaSPAdes failed!")
        logging.error(f"Return Code: {e.returncode}")
        logging.error(f"Standard Output:\n{e.stdout}")
        logging.error(f"Standard Error:\n{e.stderr}")
        raise


def validate_inputs(read1: Path, read2: Path):
    """Validate that input FASTQ files exist and are not empty."""
    for r in [read1, read2]:
        if not r.exists() or r.stat().st_size == 0:
            logging.error(f"Input FASTQ file missing or empty: {r}")
            sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Assemble metagenomic paired-end reads using metaSPAdes.")
    parser.add_argument("--read1", required=True, help="Path to R1 FASTQ file.")
    parser.add_argument("--read2", required=True, help="Path to R2 FASTQ file.")
    parser.add_argument("--output_dir", required=True, help="Directory for metaSPAdes output.")
    parser.add_argument("--threads", type=int, default=16, help="Number of CPU threads.")
    parser.add_argument("--memory", type=int, default=64, help="Max memory (GB).")
    parser.add_argument("--log_file", default="assemble_reads.log", help="Path to the log file.")
    args = parser.parse_args()

    setup_logger(Path(args.log_file))

    read1 = Path(args.read1).resolve()
    read2 = Path(args.read2).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info(f"Starting metaSPAdes assembly:")
    logging.info(f"R1={read1}\nR2={read2}\nOutput={output_dir}")

    validate_inputs(read1, read2)

    if not check_existing_outputs(output_dir):
        run_metaspades(read1, read2, output_dir, args.threads, args.memory)

    logging.info("Pipeline finished successfully.")


if __name__ == "__main__":
    main()
