import os
import csv
import json
import argparse
from pathlib import Path
from collections import defaultdict

import torch
from Bio import SeqIO
from itertools import product

# --- Settings ---
LABELS = {
    "bacteria": 0,
    "archaea": 1,
    "virus": 2,
    "plasmid": 3,
    "fungi": 4,
    "protozoa": 5,
}

K = 4  # k-mer size
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
    if total > 0:
        vec = [v / total for v in vec]
    return vec


def compute_ambiguity_rate(seq):
    return sum(1 for c in seq if c in AMBIGUOUS_BASES) / len(seq)


def process_fasta(fasta_path, gc_content, avg_length, label):
    features = []
    contig_names = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        if len(seq) < K:
            continue
        kmers = compute_kmer_freq(seq, k=K)
        amb_rate = compute_ambiguity_rate(seq)
        log_len = torch.log1p(torch.tensor(len(seq), dtype=torch.float)).item() / 20
        full_feat = [gc_content / 100.0, avg_length / 1e6, log_len, amb_rate] + kmers
        features.append(full_feat)
        contig_names.append(record.id)
    labels = [label] * len(features)
    return contig_names, features, labels


def load_preprocessing_report(report_path):
    with open(report_path) as f:
        reader = csv.DictReader(f)
        rows = [r for r in reader if r["Status"] == "KEPT"]
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset_root", type=str, required=True, help="Path to 4CAC_dataset/")
    parser.add_argument("--log_root", type=str, required=True, help="Path to logs/preprocessing-logs/")
    parser.add_argument("--output_dir", type=str, default="processed_data", help="Where to save features.pt, labels.pt, splits.pt")
    args = parser.parse_args()

    dataset_root = Path(args.dataset_root)
    log_root = Path(args.log_root)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_features, all_labels = [], []
    contig_to_idx = {}
    contig_to_split = {}
    split_masks = defaultdict(list)
    splits = ["train", "val", "test"]
    idx_counter = 0

    for split in splits:
        for category, label in LABELS.items():
            report_path = log_root / split / f"{category}_preprocessing_report.csv"
            if not report_path.exists():
                continue
            print(f"[INFO] Processing {report_path}")
            try:
                rows = load_preprocessing_report(report_path)
            except Exception as e:
                print(f"[WARNING] Skipped {report_path}: {e}")
                continue

            cleaned_dir = dataset_root / split / category / "cleaned"
            for row in rows:
                fasta_name = row["Filename"].replace(".fna", ".cleaned.fna")
                fasta_path = cleaned_dir / fasta_name
                if not fasta_path.exists():
                    continue
                gc = float(row["GC Content (%)"])
                avg_len = float(row["Avg Length"])
                contigs, feats, labels = process_fasta(fasta_path, gc, avg_len, label)
                for contig in contigs:
                    contig_to_idx[contig] = idx_counter
                    contig_to_split[contig] = split
                    split_masks[split].append(idx_counter)
                    idx_counter += 1
                all_features.extend(feats)
                all_labels.extend(labels)

    features_tensor = torch.tensor(all_features, dtype=torch.float)
    labels_tensor = torch.tensor(all_labels, dtype=torch.long)
    num_nodes = len(all_features)
    masks = {}
    for split in splits:
        mask = torch.zeros(num_nodes, dtype=torch.bool)
        mask[split_masks[split]] = True
        masks[f"{split}_mask"] = mask

    torch.save(features_tensor, output_dir / "features.pt")
    torch.save(labels_tensor, output_dir / "labels.pt")
    torch.save(masks, output_dir / "splits.pt")
    with open(output_dir / "contig_ids.json", "w") as f:
        json.dump(contig_to_idx, f)
    with open(output_dir / "contig_split_map.json", "w") as f:
        json.dump(contig_to_split, f)

    print(f"[DONE] Saved features to {output_dir}")
    print(f"[SUMMARY] Total contigs: {features_tensor.shape[0]}")
    print(f"Feature dimension: {features_tensor.shape[1]}")
    print(f"Label distribution: {dict((i, torch.sum(labels_tensor==i).item()) for i in torch.unique(labels_tensor))}")


if __name__ == "__main__":
    main()
