import os
import csv
import json
import argparse
from pathlib import Path
from collections import defaultdict
from itertools import product

import torch
from Bio import SeqIO

LABELS = {
    "bacteria": 0,
    "archaea": 1,
    "virus": 2,
    "plasmid": 3,
    "fungi": 4,
    "protozoa": 4,
}

K = 4
KMER_LIST = [''.join(p) for p in product("ACGT", repeat=K)]
KMER_INDEX = {kmer: i for i, kmer in enumerate(KMER_LIST)}
AMBIGUOUS_BASES = set("NRYSWKMBVDH")

def compute_kmer_freq(seq, k=4):
    seq = seq.upper()
    total = len(seq) - k + 1
    vec = [0] * len(KMER_LIST)
    for i in range(total):
        kmer = seq[i:i + k]
        if all(c in "ACGT" for c in kmer):
            vec[KMER_INDEX[kmer]] += 1
    return [v / total for v in vec] if total > 0 else vec

def compute_ambiguity_rate(seq):
    return sum(1 for c in seq if c in AMBIGUOUS_BASES) / len(seq)

def process_fasta(fasta_path, gc_content, avg_length, label):
    features, labels, contig_names = [], [], []
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        if len(seq) < K:
            continue
        kmers = compute_kmer_freq(seq, k=K)
        amb_rate = compute_ambiguity_rate(seq)
        log_len = torch.log1p(torch.tensor(len(seq), dtype=torch.float)).item() / 20
        full_feat = [gc_content / 100.0, avg_length / 1e6, log_len, amb_rate] + kmers
        features.append(full_feat)
        labels.append(label)
        contig_names.append(record.id)
    return contig_names, features, labels

def load_preprocessing_report(report_path):
    with open(report_path) as f:
        return [r for r in csv.DictReader(f) if r["Status"] == "KEPT"]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset_root", required=True, type=str)
    parser.add_argument("--log_root", required=True, type=str)
    parser.add_argument("--output_dir", default="processed_chunks", type=str)
    parser.add_argument("--split", type=str, help="Data split to process (e.g. train, val, test)")
    parser.add_argument("--category", type=str, help="Genome category (e.g. bacteria, archaea)")
    args = parser.parse_args()

    dataset_root = Path(args.dataset_root)
    log_root = Path(args.log_root)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    splits = [args.split] if args.split else ["train", "val", "test"]
    categories = [args.category] if args.category else list(LABELS.keys())

    for split in splits:
        for category in categories:
            label = LABELS.get(category)
            if label is None:
                print(f"[WARNING] Unknown category '{category}' — skipping.")
                continue

            print(f"[INFO] Processing {split}/{category}")
            report_path = log_root / split / f"{category}_preprocessing_report.csv"
            if not report_path.exists():
                print(f"[WARNING] Report not found: {report_path}")
                continue

            try:
                rows = load_preprocessing_report(report_path)
            except Exception as e:
                print(f"[WARNING] Failed to read {report_path}: {e}")
                continue

            cleaned_dir = dataset_root / split / category / "cleaned"
            chunk_features, chunk_labels, chunk_contigs = [], [], []
            for row in rows:
                fasta_name = row["Filename"].replace(".fna", ".cleaned.fna")
                fasta_path = cleaned_dir / fasta_name
                if not fasta_path.exists():
                    continue
                gc = float(row["GC Content (%)"])
                avg_len = float(row["Avg Length"])
                contigs, feats, labels = process_fasta(fasta_path, gc, avg_len, label)
                chunk_features.extend(feats)
                chunk_labels.extend(labels)
                chunk_contigs.extend(contigs)

            if not chunk_features:
                print(f"[WARNING] No features extracted for {split}/{category}")
                continue

            torch.save(torch.tensor(chunk_features, dtype=torch.float), output_dir / f"{split}_{category}_features.pt")
            torch.save(torch.tensor(chunk_labels, dtype=torch.long), output_dir / f"{split}_{category}_labels.pt")
            with open(output_dir / f"{split}_{category}_contigs.json", "w") as f:
                json.dump(chunk_contigs, f)

            print(f"[DONE] Saved {len(chunk_features)} contigs for {split}/{category}")

if __name__ == "__main__":
    main()