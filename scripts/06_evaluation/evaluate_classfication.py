#!/usr/bin/env python3
"""
Evaluation Pipeline for Contig Classification
---------------------------------------------
Evaluates contig classification performance (4CAC or custom system):
1. Aligns contigs to reference genomes using minimap2.
2. Assigns true classes based on alignment (≥80% identity, ≥80% coverage).
3. Compares predictions to true labels, computing precision, recall, and F1 per class.
4. Generates a confusion matrix in CSV format.

Features:
- Logs to both console and file.
- Validates input files before execution.
- Clear, reproducible outputs suitable for publications.
"""

import argparse
import subprocess
import logging
from pathlib import Path
from collections import defaultdict
from typing import Dict
import csv
import sys
from scripts.utils.logger import setup_logger


def validate_inputs(*files: Path):
    """Ensure all required input files exist and are non-empty."""
    for f in files:
        if not f.exists() or f.stat().st_size == 0:
            logging.error(f"Required input file missing or empty: {f}")
            sys.exit(1)


def run_minimap2(reference: Path, contigs: Path, output_paf: Path, threads: int = 8):
    """Run minimap2 to align contigs to reference genomes."""
    cmd = ["minimap2", "-x", "asm5", "-t", str(threads), str(reference), str(contigs)]
    logging.info(f"Running minimap2: {' '.join(cmd)}")
    try:
        with open(output_paf, "w") as paf_out:
            subprocess.run(cmd, stdout=paf_out, check=True)
        logging.info(f"Minimap2 alignment complete. Output: {output_paf}")
    except subprocess.CalledProcessError as e:
        logging.error("Minimap2 failed!")
        sys.exit(1)


def parse_paf(paf_file: Path, ref_class_map: Dict[str, str],
              min_identity: float = 0.8, min_coverage: float = 0.8) -> Dict[str, str]:
    """Parse PAF file to assign true classes to contigs."""
    contig_class = {}
    for line in open(paf_file):
        fields = line.strip().split("\t")
        if len(fields) < 12:
            continue
        query, qlen, qstart, qend, _, target, tlen, tstart, tend, match, block, _ = fields[:12]
        qlen, qstart, qend, match, block = map(int, (qlen, qstart, qend, match, block))

        identity = match / block if block > 0 else 0
        coverage = (qend - qstart) / qlen if qlen > 0 else 0

        if identity >= min_identity and coverage >= min_coverage:
            true_class = ref_class_map.get(target)
            if true_class:
                if query in contig_class and contig_class[query] != true_class:
                    contig_class[query] = "ambiguous"
                else:
                    contig_class[query] = true_class

    filtered = {c: cl for c, cl in contig_class.items() if cl != "ambiguous"}
    logging.info(f"True class labels assigned to {len(filtered)} contigs (ambiguous excluded).")
    return filtered


def load_reference_class_map(reference_fasta: Path) -> Dict[str, str]:
    """Extract class labels from reference fasta headers (>contig_id|class)."""
    ref_map = {}
    for line in open(reference_fasta):
        if line.startswith(">"):
            parts = line[1:].strip().split("|")
            ref_map[parts[0]] = parts[-1]
    logging.info(f"Loaded class labels for {len(ref_map)} reference sequences.")
    return ref_map


def load_predicted_classes(pred_file: Path) -> Dict[str, str]:
    """Load predicted classes (4CAC fasta or generic TSV)."""
    pred_map = {}
    if pred_file.suffix in [".fasta", ".fa"]:
        for line in open(pred_file):
            if line.startswith(">"):
                parts = line[1:].strip().split()
                pred_map[parts[0]] = parts[1] if len(parts) > 1 else "unclassified"
    else:
        for line in open(pred_file):
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                pred_map[parts[0]] = parts[1]
    logging.info(f"Loaded predictions for {len(pred_map)} contigs.")
    return pred_map


def compute_metrics(true_classes: Dict[str, str], pred_classes: Dict[str, str], confusion_file: Path):
    """Compute precision, recall, F1 per class and write a confusion matrix."""
    classes = sorted(set(true_classes.values()) | set(pred_classes.values()))
    tp, fp, fn = defaultdict(int), defaultdict(int), defaultdict(int)
    confusion = defaultdict(lambda: defaultdict(int))

    for contig, true_c in true_classes.items():
        pred_c = pred_classes.get(contig, "unclassified")
        confusion[true_c][pred_c] += 1
        if pred_c == true_c:
            tp[true_c] += 1
        else:
            fp[pred_c] += 1
            fn[true_c] += 1

    logging.info("=== Evaluation Metrics ===")
    for cls in classes:
        precision = tp[cls] / (tp[cls] + fp[cls]) if (tp[cls] + fp[cls]) else 0
        recall = tp[cls] / (tp[cls] + fn[cls]) if (tp[cls] + fn[cls]) else 0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0
        logging.info(f"{cls}: Precision={precision:.3f}, Recall={recall:.3f}, F1={f1:.3f}")

    with open(confusion_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        header = ["True\\Pred"] + classes
        writer.writerow(header)
        for true_c in classes:
            row = [true_c] + [confusion[true_c][pred_c] for pred_c in classes]
            writer.writerow(row)
    logging.info(f"Confusion matrix saved to {confusion_file}")


def main():
    parser = argparse.ArgumentParser(description="Evaluate classification performance.")
    parser.add_argument("--reference", required=True, help="Reference fasta with class labels in headers.")
    parser.add_argument("--contigs", required=True, help="Assembled contigs fasta.")
    parser.add_argument("--predictions", required=True, help="Predicted classes (FASTA or TSV).")
    parser.add_argument("--output_paf", default="alignment.paf", help="Minimap2 alignment output.")
    parser.add_argument("--threads", type=int, default=8, help="Threads for minimap2.")
    parser.add_argument("--confusion_file", default="confusion_matrix.csv", help="Confusion matrix CSV output.")
    parser.add_argument("--log_file", default="evaluate.log", help="Log file.")
    args = parser.parse_args()

    setup_logger(Path(args.log_file))

    reference, contigs, predictions = Path(args.reference).resolve(), Path(args.contigs).resolve(), Path(args.predictions).resolve()
    output_paf, confusion_file = Path(args.output_paf).resolve(), Path(args.confusion_file).resolve()

    logging.info("=== Starting Evaluation Pipeline ===")
    logging.info(f"Reference: {reference}\nContigs: {contigs}\nPredictions: {predictions}\nThreads: {args.threads}")

    validate_inputs(reference, contigs, predictions)

    run_minimap2(reference, contigs, output_paf, args.threads)

    true_classes = parse_paf(output_paf, load_reference_class_map(reference))
    pred_classes = load_predicted_classes(predictions)

    compute_metrics(true_classes, pred_classes, confusion_file)
    logging.info("Evaluation complete.")

if __name__ == "__main__":
    main()
