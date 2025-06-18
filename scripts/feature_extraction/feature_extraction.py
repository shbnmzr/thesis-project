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
    "viral": 2,
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


def load_preprocessing_report(report_path):
    with open(report_path) as f:
        return [r for r in csv.DictReader(f) if r["Status"] == "KEPT"]


def process_fasta_batchwise(fasta_path, gc_content, avg_length, label, batch_size, batch_callback):
    batch_features, batch_labels, batch_contigs = [], [], []
    count = 0

    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        if len(seq) < K:
            continue

        kmers = compute_kmer_freq_jellyfish(seq)
        amb_rate = compute_ambiguity_rate(seq)
        log_len = torch.log1p(torch.tensor(len(seq), dtype=torch.float)).item() / 20
        full_feat = [gc_content / 100.0, avg_length / 1e6, log_len, amb_rate] + kmers

        batch_features.append(full_feat)
        batch_labels.append(label)
        batch_contigs.append(record.id)
        count += 1

        if count % batch_size == 0:
            batch_callback(batch_features, batch_labels, batch_contigs)
            batch_features, batch_labels, batch_contigs = [], [], []

    if batch_features:
        batch_callback(batch_features, batch_labels, batch_contigs)


def merge_batches(output_dir, split, category):
    feature_parts = sorted(output_dir.glob(f"{split}_{category}_features_part*.pt"))
    label_parts = sorted(output_dir.glob(f"{split}_{category}_labels_part*.pt"))
    contig_parts = sorted(output_dir.glob(f"{split}_{category}_contigs_part*.json"))

    all_features = [torch.load(p, weights_only=True) for p in feature_parts]
    all_labels = [torch.load(p, weights_only=True) for p in label_parts]
    all_contigs = []
    for p in contig_parts:
        with open(p) as f:
            all_contigs.extend(json.load(f))

    merged_features = torch.cat(all_features)
    merged_labels = torch.cat(all_labels)

    torch.save(merged_features, output_dir / f"{split}_{category}_features.pt")
    torch.save(merged_labels, output_dir / f"{split}_{category}_labels.pt")
    with open(output_dir / f"{split}_{category}_contigs.json", "w") as f:
        json.dump(all_contigs, f)

    print(f"[MERGED] {split}/{category}: {len(merged_features)} contigs")

    # Clean up part files
    for p in feature_parts + label_parts + contig_parts:
        try:
            p.unlink()
            print(f"[DELETED] {p.name}")
        except Exception as e:
            print(f"[ERROR] Failed to delete {p.name}: {e}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset_root", required=True)
    parser.add_argument("--log_root", required=True)
    parser.add_argument("--output_dir", default="processed_chunks")
    parser.add_argument("--split", type=str)
    parser.add_argument("--category", type=str)
    parser.add_argument("--batch_size", type=int, default=5000)
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
            part_id = 0

            def save_batch(features, labels, contigs):
                nonlocal part_id
                torch.save(torch.tensor(features, dtype=torch.float), output_dir / f"{split}_{category}_features_part{part_id}.pt")
                torch.save(torch.tensor(labels, dtype=torch.long), output_dir / f"{split}_{category}_labels_part{part_id}.pt")
                with open(output_dir / f"{split}_{category}_contigs_part{part_id}.json", "w") as f:
                    json.dump(contigs, f)
                print(f"[SAVED] {split}/{category} - batch {part_id} with {len(features)} contigs")
                part_id += 1

            for row in rows:
                fasta_name = row["Filename"].replace(".fna", ".cleaned.fna")
                fasta_path = cleaned_dir / fasta_name
                if not fasta_path.exists():
                    continue
                gc = float(row["GC Content (%)"])
                avg_len = float(row["Avg Length"])
                process_fasta_batchwise(fasta_path, gc, avg_len, label, args.batch_size, save_batch)

            # Merge parts
            merge_batches(output_dir, split, category)


if __name__ == "__main__":
    main()
