import os
import csv
import json
import argparse
import subprocess
import tempfile
from pathlib import Path
from itertools import product

import numpy as np
import torch
from Bio import SeqIO

# --- Constants ---
K = 4
KMER_LIST = [''.join(p) for p in product("ACGT", repeat=K)]
KMER_INDEX = {kmer: i for i, kmer in enumerate(KMER_LIST)}
AMBIGUOUS_BASES = set("NRYSWKMBVDH")
LABELS = {
    "bacteria": 0,
    "archaea": 1,
    "virus": 2,
    "plasmid": 3,
    "fungi": 4,
    "protozoa": 4,
}


# --- Functions ---
def compute_kmer_freq_jellyfish(seq, k=4):
    vec = np.zeros(4**k, dtype=np.float32)
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        fasta_path = tmpdir / "temp.fa"
        output_txt = tmpdir / "counts.txt"

        with open(fasta_path, "w") as f:
            f.write(">seq\n" + seq + "\n")

        subprocess.run(["bash", "count_kmers.sh", str(fasta_path), str(output_txt)], check=True)

        total = 0
        with open(output_txt) as f:
            for line in f:
                kmer, count = line.strip().split()
                if kmer in KMER_INDEX:
                    vec[KMER_INDEX[kmer]] = int(count)
                    total += int(count)
        if total > 0:
            vec /= total
    return vec.tolist()


def compute_ambiguity_rate(seq):
    return sum(1 for c in seq if c in AMBIGUOUS_BASES) / len(seq)


def process_fasta(fasta_path, gc_content, avg_length, label):
    features, labels, contig_names = [], [], []
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        if len(seq) < K:
            continue
        kmers = compute_kmer_freq_jellyfish(seq)
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
    parser.add_argument("--dataset_root", required=True)
    parser.add_argument("--log_root", required=True)
    parser.add_argument("--output_dir", default="processed_chunks")
    parser.add_argument("--split", type=str)
    parser.add_argument("--category", type=str)
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
                print(f"[WARNING] Unknown category '{category}'")
                continue

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
                continue

            torch.save(torch.tensor(chunk_features, dtype=torch.float), output_dir / f"{split}_{category}_features.pt")
            torch.save(torch.tensor(chunk_labels, dtype=torch.long), output_dir / f"{split}_{category}_labels.pt")
            with open(output_dir / f"{split}_{category}_contigs.json", "w") as f:
                json.dump(chunk_contigs, f)

            print(f"[DONE] {split}/{category}: {len(chunk_features)} contigs")


if __name__ == "__main__":
    main()
