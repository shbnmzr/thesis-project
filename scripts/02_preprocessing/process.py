#!/usr/bin/env python3
"""
Genome Preprocessing and Deduplication
----------------------------------------------------------------------------------

Features:
1. Cleans genomes by:
   - Minimum contig length
   - Ambiguous base filtering
   - Removes duplicate contigs within a genome

2. Deduplicates across splits (train/val/test):
   - Exact duplicates (MD5 hash)
   - Near-identical genomes (Mash distance ≤ 0.05 ≈ ANI ≥ 95%)
   - Compares all pairs: train-val, train-test, val-test
   - Batch-wise Mash for memory efficiency and long argument list prevention

Outputs:
- Cleaned FASTAs under train/val/test/<category>/cleaned
- Preprocessing and deduplication reports in the log directory
"""

import argparse
import logging
import sys
import hashlib
import subprocess
from pathlib import Path
from Bio import SeqIO
import csv

# --- Ambiguous base set ---
AMBIGUOUS_BASES = set("NRYSWKMBVDH")

# -----------------------------
# Logging
# -----------------------------
def setup_logging(log_dir: Path, category: str):
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

# -----------------------------
# Cleaning Functions
# -----------------------------
def is_ambiguous(seq):
    return sum(base in AMBIGUOUS_BASES for base in seq)

def safe_seqio_parse(fasta_path):
    try:
        return list(SeqIO.parse(fasta_path, "fasta"))
    except Exception as e:
        logging.error(f"Failed to parse {fasta_path.name}: {e}")
        return []

def clean_fasta(fasta_path, cleaned_dir, min_len, max_ambig):
    """Clean a genome FASTA by removing short/ambiguous/duplicate contigs."""
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

def get_gc_and_avg_len(fasta_path):
    seqs = safe_seqio_parse(fasta_path)
    if not seqs:
        return 0.0, 0.0
    total_gc = sum((str(rec.seq).upper().count("G") + str(rec.seq).upper().count("C")) for rec in seqs)
    total_len = sum(len(rec.seq) for rec in seqs)
    return (100 * total_gc / total_len, total_len / len(seqs))

def compute_genome_hash(fasta_path):
    """Compute a hash for the entire genome (concatenated contigs)."""
    seqs = safe_seqio_parse(fasta_path)
    concat_seq = "".join(sorted([str(rec.seq).upper() for rec in seqs]))
    return hashlib.md5(concat_seq.encode()).hexdigest()

# -----------------------------
# Deduplication Functions
# -----------------------------
def deduplicate_across_splits(base_dir: Path, category: str, log_dir: Path):
    """
    Remove genomes from val/test if identical to train genomes (MD5 hash).
    Also remove duplicates between val and test.
    """
    logging.info("Starting cross-split deduplication (exact duplicates)...")
    splits = ["train", "val", "test"]
    genome_hashes = {s: {} for s in splits}

    for split in splits:
        cleaned_dir = base_dir / split / category / "cleaned"
        if not cleaned_dir.exists():
            continue
        for f in cleaned_dir.glob("*.fna"):
            genome_hashes[split][f] = compute_genome_hash(f)

    removed = []
    train_hashes = set(genome_hashes["train"].values())
    val_hashes = set(genome_hashes["val"].values())

    for split in ["val", "test"]:
        for f, h in list(genome_hashes[split].items()):
            if h in train_hashes and f.exists():
                logging.warning(f"Removing duplicate genome in {split}: {f.name}")
                f.unlink()
                removed.append((split, f.name, f"Exact duplicate (MD5, Train vs {split.capitalize()})"))

    # Val vs Test
    for f, h in list(genome_hashes["test"].items()):
        if h in val_hashes and f.exists():
            logging.warning(f"Removing duplicate genome in test: {f.name}")
            f.unlink()
            removed.append(("test", f.name, "Exact duplicate (MD5, Val vs Test)"))

    # Save report
    report_file = log_dir / f"{category}_deduplication_report.csv"
    with open(report_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Split", "Filename_Removed", "Reason"])
        for row in removed:
            writer.writerow(row)
    logging.info(f"Exact deduplication complete. Report saved: {report_file}")

def mash_compare(query_sketch, ref_sketch):
    """Run Mash distance between query and reference sketches."""
    cmd_dist = ["mash", "dist", ref_sketch, query_sketch]
    result = subprocess.run(cmd_dist, capture_output=True, text=True, check=True)
    distances = {}
    for line in result.stdout.strip().split("\n"):
        if not line.strip():
            continue
        ref, query, dist, _, _ = line.split("\t")
        distances[(Path(query).name, Path(ref).name)] = float(dist)
    return distances

def sketch_genomes_in_batches(genome_paths, sketch_prefix, batch_size=200):
    """Sketch genomes in batches and merge into one .msh file."""
    batch_files = []
    for i in range(0, len(genome_paths), batch_size):
        batch = genome_paths[i:i + batch_size]
        batch_sketch = f"{sketch_prefix}_batch_{i//batch_size}"
        cmd_sketch = ["mash", "sketch", "-o", batch_sketch] + [str(g) for g in batch]
        logging.info(f"Sketching batch {i//batch_size + 1}: {len(batch)} genomes → {batch_sketch}.msh")
        subprocess.run(cmd_sketch, check=True)
        batch_files.append(f"{batch_sketch}.msh")

    if len(batch_files) == 1:
        return batch_files[0]

    merged_sketch = f"{sketch_prefix}_merged"
    cmd_paste = ["mash", "paste", merged_sketch] + batch_files
    logging.info(f"Merging {len(batch_files)} batch sketches → {merged_sketch}.msh")
    subprocess.run(cmd_paste, check=True)
    return f"{merged_sketch}.msh"

def remove_near_identicals_batch(base_dir: Path, category: str, log_dir: Path,
                                 threshold=0.05, batch_size=100):
    """
    Remove near-identical genomes (Mash distance ≤ threshold) across:
    - train vs val
    - train vs test
    - val vs test
    """
    logging.info(f"Running Mash-based near-identical removal (batch size={batch_size})...")

    def compare_and_remove(ref_split, query_split):
        ref_paths = list((base_dir / ref_split / category / "cleaned").glob("*.fna"))
        query_paths = list((base_dir / query_split / category / "cleaned").glob("*.fna"))
        if not ref_paths or not query_paths:
            return []

        ref_sketch = sketch_genomes_in_batches(ref_paths, str(log_dir / f"{category}_{ref_split}_mash"), batch_size=200)
        removed_local = []

        for i in range(0, len(query_paths), batch_size):
            batch = query_paths[i:i + batch_size]
            batch_sketch = str(log_dir / f"{category}_{query_split}_batch_{i//batch_size}")
            cmd_sketch = ["mash", "sketch", "-o", batch_sketch] + [str(g) for g in batch]
            subprocess.run(cmd_sketch, check=True)
            distances = mash_compare(f"{batch_sketch}.msh", ref_sketch)

            for genome in batch:
                for ref_file in ref_paths:
                    key = (genome.name, ref_file.name)
                    if key in distances and distances[key] <= threshold:
                        if genome.exists():
                            logging.warning(f"Removing near-identical genome in {query_split}: {genome.name}")
                            genome.unlink()
                            removed_local.append((
                                query_split,
                                genome.name,
                                f"Mash_dist={distances[key]:.4f} ({ref_split.capitalize()} vs {query_split.capitalize()})"
                            ))
                        break
        return removed_local

    removed = []
    removed += compare_and_remove("train", "val")
    removed += compare_and_remove("train", "test")
    removed += compare_and_remove("val", "test")

    report_file = log_dir / f"{category}_near_identical_report.csv"
    with open(report_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Split", "Filename_Removed", "Reason"])
        for row in removed:
            writer.writerow(row)
    logging.info(f"Near-identical removal complete. Report saved: {report_file}")

# -----------------------------
# Main
# -----------------------------
def main():
    parser = argparse.ArgumentParser(description="Clean and deduplicate FASTA files by length, ambiguity, and similarity.")
    parser.add_argument("--base_dir", type=str, required=True, help="Base directory of the dataset")
    parser.add_argument("--category", type=str, required=True, help="Genome category (e.g., archaea)")
    parser.add_argument("--min_len", type=int, default=1000, help="Minimum contig length")
    parser.add_argument("--max_ambig", type=float, default=0.05, help="Maximum allowed ambiguous base ratio")
    parser.add_argument("--min_contigs", type=int, default=10, help="Minimum number of contigs to keep genome")
    parser.add_argument("--log_dir", type=str, required=True, help="Directory to store logs")
    args = parser.parse_args()

    base_path = Path(args.base_dir)
    splits = ["train", "val", "test"]
    log_dir = Path(args.log_dir)

    setup_logging(log_dir, args.category)

    for split in splits:
        genome_dir = base_path / split / args.category
        cleaned_dir = genome_dir / "cleaned"

        if not genome_dir.exists():
            logging.warning(f"Genome directory not found: {genome_dir.resolve()}")
            continue

        cleaned_dir.mkdir(parents=True, exist_ok=True)
        fasta_files = [f for ext in ["*.fna", "*.fa", "*.fasta"] for f in genome_dir.glob(ext)]
        if not fasta_files:
            logging.warning(f"No FASTA files found in: {genome_dir.resolve()}")
            continue

        logging.info(f"Starting preprocessing for category '{args.category}' [{split}]...")
        report = []

        for fasta_file in fasta_files:
            logging.info(f"Processing {fasta_file.name}")
            cleaned_path, before, after = clean_fasta(fasta_file, cleaned_dir, args.min_len, args.max_ambig)

            if cleaned_path and after >= args.min_contigs:
                gc, avg_len = get_gc_and_avg_len(cleaned_path)
                status = "KEPT"
            elif cleaned_path:
                if cleaned_path.exists():
                    cleaned_path.unlink()
                gc, avg_len = 0.0, 0.0
                status = "DELETED (too few contigs)"
            else:
                gc, avg_len = 0.0, 0.0
                status = "DELETED (empty or error)"

            report.append((fasta_file.name, before, after, f"{gc:.2f}", f"{avg_len:.1f}", status))

        report_file = log_dir / f"{args.category}_{split}_preprocessing_report.csv"
        with open(report_file, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Filename", "Original Contigs", "After Filtering", "GC Content (%)", "Avg Length", "Status"])
            for row in report:
                writer.writerow(row)
        logging.info(f"Preprocessing complete for {split}. Report: {report_file}")

    deduplicate_across_splits(base_path, args.category, log_dir)
    remove_near_identicals_batch(base_path, args.category, log_dir, threshold=0.05, batch_size=100)

if __name__ == "__main__":
    main()
