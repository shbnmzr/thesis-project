import os
import requests
import gzip
import shutil
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO

CATEGORY = "viral"  # Change as needed
BASE_DIR = Path("../4CAC_dataset")
METADATA_DIR = BASE_DIR / "metadata"
TRAIN_DIR = BASE_DIR / "train" / CATEGORY
VAL_DIR = BASE_DIR / "val" / CATEGORY
TEST_DIR = BASE_DIR / "test" / CATEGORY

TRAIN_CUTOFF = datetime.strptime("2021-12-01", "%Y-%m-%d")
VAL_CUTOFF = datetime.strptime("2022-04-01", "%Y-%m-%d")
TEST_END = datetime.strptime("2023-04-01", "%Y-%m-%d")
MAX_WORKERS = 8

for d in [TRAIN_DIR, VAL_DIR, TEST_DIR, METADATA_DIR]:
    d.mkdir(parents=True, exist_ok=True)

assembly_summary_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/genbank/{CATEGORY}/assembly_summary.txt"
summary_file = METADATA_DIR / f"{CATEGORY}_assembly_summary.txt"

print(f"[INFO] Downloading GenBank assembly summary for '{CATEGORY}'...")
response = requests.get(assembly_summary_url)
response.raise_for_status()
with open(summary_file, "w") as f:
    f.write(response.text)

with open(summary_file) as f:
    header = next(line for line in f if line.startswith("#assembly_accession")).strip().split('\t')
ftp_col_index = header.index("ftp_path")
date_col_index = header.index("seq_rel_date")

dataset_entries = {"train": [], "val": [], "test": []}

with open(summary_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split('\t')
        if len(parts) > 14:
            try:
                release_date = datetime.strptime(parts[14], "%Y-%m-%d")
                if parts[11] in {"Complete Genome", "Chromosome"} and parts[10] == "latest":
                    path = parts[ftp_col_index]
                    if release_date < TRAIN_CUTOFF:
                        dataset_entries["train"].append((path, release_date))
                    elif release_date < VAL_CUTOFF:
                        dataset_entries["val"].append((path, release_date))
                    elif release_date < TEST_END:
                        dataset_entries["test"].append((path, release_date))
            except ValueError:
                continue

for split in ["train", "val", "test"]:
    print(f"[INFO] {split.upper()}: {len(dataset_entries[split])} genomes")
    with open(METADATA_DIR / f"{CATEGORY}_{split}_ftp_paths.txt", "w") as f:
        for path, _ in dataset_entries[split]:
            f.write(path + "\n")

# FASTA validation helpers
def is_valid_fasta(path):
    try:
        with open(path) as f:
            for line in f:
                if line.startswith(">"):
                    return True
        return False
    except:
        return False

def has_valid_contig(path, min_len=1000):
    try:
        for rec in SeqIO.parse(path, "fasta"):
            if len(rec.seq) >= min_len:
                return True
        return False
    except:
        return False

def download_file(url, dest_path):
    try:
        with requests.get(url.replace("ftp://", "https://"), stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(dest_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        return True
    except Exception:
        return False

def decompress_gz(gz_file):
    try:
        output_path = gz_file.with_suffix('')
        with gzip.open(gz_file, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        gz_file.unlink()
        if output_path.stat().st_size < 500 or not is_valid_fasta(output_path) or not has_valid_contig(output_path):
            print(f"[WARNING] Removing invalid or tiny FASTA: {output_path.name}")
            output_path.unlink()
            return False
        return True
    except Exception as e:
        print(f"[WARNING] Failed to decompress: {gz_file.name} - {e}")
        return False

def download_split(split_name, entries, output_dir):
    print(f"[INFO] Downloading {split_name.upper()} genomes...")
    failed_urls = []
    output_dir.mkdir(parents=True, exist_ok=True)

    def download_worker(path):
        accession = path.split("/")[-1]
        url = f"{path}/{accession}_genomic.fna.gz"
        filename = url.split("/")[-1]
        dest_path = output_dir / filename
        if dest_path.exists():
            return 'skipped'
        if download_file(url, dest_path):
            return 'downloaded'
        return url  # failed

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(download_worker, path): path for path, _ in entries}
        for future in as_completed(futures):
            result = future.result()
            if result not in {"downloaded", "skipped"}:
                failed_urls.append(result)

    if failed_urls:
        with open(METADATA_DIR / f"{CATEGORY}_{split_name}_failed.txt", "w") as f:
            for url in failed_urls:
                f.write(url + "\n")

    print(f"[INFO] Decompressing {split_name.upper()} files...")
    gz_files = list(output_dir.glob("*.gz"))
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        list(executor.map(decompress_gz, gz_files))

for split, entries in dataset_entries.items():
    split_dir = BASE_DIR / split / CATEGORY
    download_split(split, entries, split_dir)

print(f"[INFO] Completed downloading and organizing datasets for category: {CATEGORY}")