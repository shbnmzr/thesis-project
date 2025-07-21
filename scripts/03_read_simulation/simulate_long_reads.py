#!/usr/bin/env python3
"""
Simulate Nanopore long reads

Features:
1. Uses NanoSim to generate realistic ONT-like long reads.
2. Supports multiple scenarios (generic + filtered).
3. Outputs FASTQ files suitable for metagenomic assembly (metaFlye).
4. Saves metadata for reproducibility.

Example:
    python simulate_long_reads_4cac.py \
        --fasta ./4CAC_dataset/test_prep/all_references.fasta \
        --abundance_generic ./4CAC_dataset/test_prep/abundance_generic.txt \
        --abundance_filtered ./4CAC_dataset/test_prep/abundance_filtered.txt \
        --output_dir ./4CAC_dataset/simulated_long \
        --n_reads 500000 \
        --model_path ./nanosim_models/ecoli_R9_1D
"""

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path
from scripts.utils.logger import setup_logger  # your custom logger

# -----------------------------
# NANO SIM WRAPPER
# -----------------------------
def run_nanosim(fasta: Path, abundance: Path, output_prefix: Path,
                n_reads: int, model_path: Path):
    """
    Run NanoSim to simulate Nanopore long reads.

    Args:
        fasta (Path): Combined reference FASTA file.
        abundance (Path): Abundance file for this scenario.
        output_prefix (Path): Prefix for output FASTQs.
        n_reads (int): Total number of reads to simulate.
        model_path (Path): Path to NanoSim pretrained model directory.
    """
    cmd = [
        "nanosim-h",
        "-DNA", str(fasta),
        "--abundance", str(abundance),
        "--model_path", str(model_path),
        "-n", str(n_reads),
        "--output", str(output_prefix)
    ]
    logging.info(f"Running NanoSim: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
        logging.info(f"Simulation complete: {output_prefix}.fastq")
    except subprocess.CalledProcessError as e:
        logging.error(f"NanoSim failed: {e}")
        sys.exit(1)
    except FileNotFoundError:
        logging.error("NanoSim is not installed or not in PATH.")
        sys.exit(1)

# -----------------------------
# MAIN
# -----------------------------
def main():
    parser = argparse.ArgumentParser(description="Simulate Nanopore long reads for metagenomic test sets.")
    parser.add_argument("--fasta", required=True, help="Combined reference FASTA (with class prefixes).")
    parser.add_argument("--abundance_generic", required=True, help="Abundance file for generic scenario.")
    parser.add_argument("--abundance_filtered", required=True, help="Abundance file for filtered scenario.")
    parser.add_argument("--output_dir", required=True, help="Directory to save simulated long reads.")
    parser.add_argument("--model_path", required=True, help="NanoSim pretrained model directory.")
    parser.add_argument("--n_reads", type=int, default=500000, help="Total reads per scenario.")
    parser.add_argument("--log_file", default="./logs/simulate_long_reads.log",
                        help="Log file (default: ./logs/simulate_long_reads.log)")
    args = parser.parse_args()

    # Setup logging
    setup_logger(Path(args.log_file))

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Run generic simulation
    run_nanosim(Path(args.fasta), Path(args.abundance_generic),
                output_dir / "generic_nanopore", args.n_reads, Path(args.model_path))

    # Run filtered simulation
    run_nanosim(Path(args.fasta), Path(args.abundance_filtered),
                output_dir / "filtered_nanopore", args.n_reads, Path(args.model_path))

    # Save metadata
    metadata = {
        "fasta": str(args.fasta),
        "abundance_generic": str(args.abundance_generic),
        "abundance_filtered": str(args.abundance_filtered),
        "model_path": str(args.model_path),
        "n_reads": args.n_reads
    }
    with open(output_dir / "simulation_metadata.json", "w") as f:
        json.dump(metadata, f, indent=4)

    logging.info(f"Metadata saved: {output_dir / 'simulation_metadata.json'}")
    logging.info("Long-read simulation complete for both scenarios.")


if __name__ == "__main__":
    main()
