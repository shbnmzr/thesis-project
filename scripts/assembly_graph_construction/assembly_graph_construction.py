#!/usr/bin/env python3

import subprocess
from pathlib import Path
import argparse
import logging
import shutil
import time
import os

# === Configuration ===
ART_BIN = "art_illumina"
SPADES_BIN = "/path/to/spades/bin/spades.py"
FASTA_SUFFIX = ".cleaned.fna"
DEFAULT_COVERAGE = "10"
TMP_BASE = Path(os.environ.get("TMPDIR", "/tmp"))


def setup_logging(log_file: Path = None):
    """Configure logging output."""
    handlers = [logging.StreamHandler()]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
    )


def simulate_reads(fna_path: Path, output_prefix: Path, coverage: str):
    """Simulate paired-end reads from a FASTA using ART with lower coverage."""
    cmd = [
        ART_BIN, "-ss", "HS25", "-i", str(fna_path), "-p",
        "-l", "150", "-f", coverage, "-m", "200", "-s", "10",
        "-o", str(output_prefix)
    ]
    subprocess.run(cmd, check=True)


def run_spades(read1: Path, read2: Path, output_dir: Path, tmp_dir: Path):
    """Run SPAdes using meta mode with temporary output location."""
    cmd = [
        SPADES_BIN,
        "-1", str(read1), "-2", str(read2),
        "-o", str(output_dir),
        "--meta",
        "--only-assembler",
        "--tmp-dir", str(tmp_dir)
    ]
    subprocess.run(cmd, check=True)


def cleanup_intermediates(output_dir: Path):
    """Remove everything except the final .gfa graph."""
    for item in output_dir.iterdir():
        if item.name != "assembly_graph_with_scaffolds.gfa":
            if item.is_dir():
                shutil.rmtree(item, ignore_errors=True)
            else:
                item.unlink(missing_ok=True)


def process_genome(fna_path: Path, output_dir: Path, coverage: str):
    """Process a single genome: simulate reads, run SPAdes, and clean."""
    genome_id = fna_path.stem.replace(".cleaned", "")
    genome_out_dir = output_dir / genome_id
    genome_out_dir.mkdir(parents=True, exist_ok=True)

    gfa_path = genome_out_dir / "assembly_graph_with_scaffolds.gfa"
    if gfa_path.exists():
        logging.info(f"[SKIP] GFA exists for {genome_id}")
        return

    read_prefix = genome_out_dir / genome_id
    read1 = read_prefix.with_name(read_prefix.name + "1.fq")
    read2 = read_prefix.with_name(read_prefix.name + "2.fq")
    tmp_dir = TMP_BASE / f"tmp_spades_{genome_id}"

    try:
        tmp_dir.mkdir(parents=True, exist_ok=True)

        logging.info(f"[SIMULATE] {genome_id}")
        simulate_reads(fna_path, read_prefix, coverage)

        logging.info(f"[SPADES] {genome_id}")
        start = time.time()
        run_spades(read1, read2, genome_out_dir, tmp_dir)
        elapsed = time.time() - start
        logging.info(f"[DONE] {genome_id} in {elapsed:.1f} sec")

        read1.unlink(missing_ok=True)
        read2.unlink(missing_ok=True)
        cleanup_intermediates(genome_out_dir)

    except subprocess.CalledProcessError as e:
        logging.error(f"[ERROR] SPAdes or ART failed for {genome_id}: {e}")
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def main(split: str, category: str, input_root: Path, output_root: Path, coverage: str):
    """Main orchestration logic."""
    input_dir = input_root / split / category / "cleaned"
    output_dir = output_root / split / category
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_dir.exists():
        logging.error(f"[ABORT] Input directory not found: {input_dir}")
        return

    fna_files = sorted(input_dir.glob(f"*{FASTA_SUFFIX}"))
    logging.info(f"Found {len(fna_files)} genomes in '{category}' split '{split}'")

    for fna_path in fna_files:
        process_genome(fna_path, output_dir, coverage)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Efficiently generate SPAdes assembly graphs.")
    parser.add_argument("--split", type=str, required=True, help="Dataset split: train / val / test")
    parser.add_argument("--category", type=str, required=True, help="Taxonomic group name (e.g., archaea)")
    parser.add_argument("--input_root", type=Path, default=Path("/nfsd/bcb/bcbg/zareshahraki/4CAC_dataset"),
                        help="Root input directory containing preprocessed contigs")
    parser.add_argument("--output_root", type=Path, default=Path("/nfsd/bcb/bcbg/zareshahraki/assembly_graphs"),
                        help="Directory to store SPAdes graphs")
    parser.add_argument("--coverage", type=str, default=DEFAULT_COVERAGE,
                        help="Fold coverage for ART read simulation (e.g., 5)")
    args = parser.parse_args()

    setup_logging()
    main(args.split, args.category, args.input_root, args.output_root, args.coverage)
