#!/usr/bin/env python3
"""
Preprocess genome FASTA files by filtering contigs based on:
- Minimum length
- Maximum ambiguous base ratio
- Minimum number of contigs per genome

Outputs:
- Cleaned FASTA files (only genomes that pass filters)
- CSV report summarizing preprocessing
- Logging of all actions
"""

import argparse
import logging
import sys
import hashlib
from pathlib import Path
from Bio import SeqIO
import csv

# -----------------------------
# GLOBALS
# -----------------------------
AMBIGUOUS_BASES = set("NRYSWKMBVDH")  # IUPAC ambiguity codes


# -----------------------------
# LOGGING
# -----------------------------
def setup_logging(log_dir: Path, category: str):
    """Configure logging to both a file and stdout."""
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{category}_preprocessing.log"
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logging.info(f"Logging initialized. Log file: {log_file.resolve()}")


# -----------------------------
# HELPER FUNCTIONS
# -----------------------------
def is_ambiguous(seq: str) -> int:
    """Count ambiguous bases in a sequence."""
    return sum(base in AMBIGUOUS_BASES for base in seq)


def safe_seqio_parse(fasta_path: Path):
    """Safely parse FASTA files, returning a list of records."""
    try:
        return list(SeqIO.parse(fasta_path, "fasta"))
    except Exception as e:
        logging.error(f"Failed to parse {fasta_path.name}: {e}")
        return []


def clean_fasta(fasta_path: Path, cleaned_dir: Path, min_len: int, max_ambig: float):
    """
    Filter contigs in a FASTA by minimum length and ambiguity.

    Returns:
        (cleaned_path, original_count, kept_count)
    """
    seen = set()
    records = []
    total_before = 0
    seqs = safe_seqio_parse(fasta_path)

    for rec in seqs:
        total_before += 1
        seq = str(rec.seq).upper()
        if not seq:
            continue
        seq_hash = hashlib.md5(seq.encode()).hexdigest()
        if seq_hash in seen:
            continue
        seen.add(seq_hash)
        if len(seq) < min_len:
            continue
        if (is_ambiguous(seq) / len(seq)) > max_ambig:
            continue
        records.append(rec)

    if records:
        cleaned_filename = fasta_path.stem + ".cleaned.fna"
        cleaned_path = cleaned_dir / cleaned_filename
        SeqIO.write(records, cleaned_path, "fasta")
        return cleaned_path, total_before, len(records)
    else:
        return None, total_before, 0


def get_gc_and_avg_len(fasta_path: Path):
    """Calculate GC content (%) and average contig length."""
    seqs = safe_seqio_parse(fasta_path)
    if not seqs:
        return 0.0, 0.0
    total_gc = sum((str(rec.seq).upper().count("G") + str(rec.seq).upper().count("C")) for rec in seqs)
    total_len = sum(len(rec.seq) for rec in seqs)
    return (100 * total_gc / total_len, total_len / len(seqs))


# -----------------------------
# MAIN FUNCTION
# -----------------------------
def main():
    parser = argparse.ArgumentParser(description="Preprocess genome FASTA files by filtering contigs.")
    parser.add_argument("--base_dir", type=str, required=True, help="Base directory of the dataset.")
    parser.add_argument("--category", type=str, required=True, help="Genome category (e.g., archaea).")
    parser.add_argument("--min_len", type=int, default=1000, help="Minimum contig length.")
    parser.add_argument("--max_ambig", type=float, default=0.05, help="Maximum allowed ambiguous base ratio.")
    parser.add_argument("--min_contigs", type=int, default=10, help="Minimum contigs required to keep genome.")
    parser.add_argument("--log_dir", type=str, required=True, help="Directory to store logs and reports.")
    args = parser.parse_args()

    base_path = Path(args.base_dir)
    genome_dir = base_path / "train" / args.category
    cleaned_dir = genome_dir / "cleaned"
    log_dir = Path(args.log_dir)

    setup_logging(log_dir, args.category)

    if not genome_dir.exists():
        logging.error(f"Genome directory not found: {genome_dir.resolve()}")
        sys.exit(1)

    cleaned_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = [f for ext in ["*.fna", "*.fa", "*.fasta"] for f in genome_dir.glob(ext)]
    if not fasta_files:
        logging.error(f"No FASTA files found in: {genome_dir.resolve()}")
        sys.exit(1)

    logging.info(f"Starting preprocessing for category '{args.category}'...")
    report = []

    for fasta_file in fasta_files:
        logging.info(f"Processing {fasta_file.name}")
        cleaned_path, before, after = clean_fasta(fasta_file, cleaned_dir, args.min_len, args.max_ambig)

        if cleaned_path and after >= args.min_contigs:
            gc, avg_len = get_gc_and_avg_len(cleaned_path)
            status = "KEPT"
        elif cleaned_path:
            cleaned_path.unlink()
            gc, avg_len = 0.0, 0.0
            status = "DELETED (too few contigs)"
        else:
            gc, avg_len = 0.0, 0.0
            status = "DELETED (empty or error)"

        report.append((fasta_file.name, before, after, f"{gc:.2f}", f"{avg_len:.1f}", status))

    # Save CSV report
    report_file = log_dir / f"{args.category}_preprocessing_report.csv"
    with open(report_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Filename", "Original Contigs", "After Filtering", "GC Content (%)", "Avg Length", "Status"])
        writer.writerows(report)

    logging.info(f"Preprocessing complete. Report saved to: {report_file.resolve()}")


if __name__ == "__main__":
    main()
