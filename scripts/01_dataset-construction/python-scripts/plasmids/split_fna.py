#!/usr/bin/env python3
"""
Split a multi-FASTA file into individual FASTA files by accession.

Usage:
    python split_fasta_by_accession.py \
        --input_fasta ../4CAC_dataset/genomes/plasmids/all_plasmids.fasta \
        --output_dir ../4CAC_dataset/train/plasmid/raw
"""

import argparse
import logging
from pathlib import Path
from Bio import SeqIO
from collections import defaultdict

# -----------------------------
# LOGGING CONFIGURATION
# -----------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()]
)


# -----------------------------
# FUNCTION DEFINITIONS
# -----------------------------
def split_fasta(input_fasta: Path, output_dir: Path) -> None:
    """
    Split a multi-FASTA into individual FASTA files grouped by accession.

    Args:
        input_fasta (Path): Path to the input multi-FASTA file.
        output_dir (Path): Directory to save individual FASTA files.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    genomes = defaultdict(list)

    # Parse sequences and group by accession
    for record in SeqIO.parse(input_fasta, "fasta"):
        accession = record.id.split("|")[0]  # Example: "GCF_XXXXX"
        genomes[accession].append(record)

    logging.info(f"Found {len(genomes)} unique accessions in {input_fasta.name}.")

    # Write each accession to its own FASTA
    for accession, records in genomes.items():
        output_path = output_dir / f"{accession}.fna"
        SeqIO.write(records, output_path, "fasta")

    logging.info(f"Saved {len(genomes)} FASTA files to {output_dir}")


# -----------------------------
# MAIN EXECUTION
# -----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split multi-FASTA into individual FASTA files by accession.")
    parser.add_argument("--input_fasta", required=True, help="Path to the input multi-FASTA file.")
    parser.add_argument("--output_dir", required=True, help="Directory to save split FASTA files.")
    args = parser.parse_args()

    split_fasta(Path(args.input_fasta), Path(args.output_dir))
