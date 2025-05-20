import os
import requests
import time
import gzip
import shutil
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

CATEGORY = "viral"
BASE_DIR = Path("../4CAC_dataset")
METADATA_DIR = BASE_DIR / "metadata"
GENOME_DIR = BASE_DIR / "genomes" / CATEGORY

# Parameters
MIN_DATE = datetime.strptime("2023-04-01", "%Y-%m-%d")
MAX_WORKERS = 8

# Create necessary directories
GENOME_DIR.mkdir(parents=True, exist_ok=True)
METADATA_DIR.mkdir(parents=True, exist_ok=True)

# Step 1: Download the assembly summary
assembly_summary_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/genbank/{CATEGORY}/assembly_summary.txt"
summary_file = METADATA_DIR / f"{CATEGORY}_assembly_summary.txt"

print(f"[INFO] Downloading assembly summary for '{CATEGORY}'...")
response = requests.get(assembly_summary_url)
response.raise_for_status()
with open(summary_file, "w") as f:
    f.write(response.text)

# Find ftp_path column index
with open(summary_file) as f:
    header = next(line for line in f if line.startswith("#assembly_accession")).strip().split('\t')
ftp_col_index = header.index("ftp_path")

# Extract relevant ftp paths (filtered by date)
ftp_paths = []
with open(summary_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split('\t')
        if len(parts) > 14:
            try:
                release_date = datetime.strptime(parts[14], "%Y-%m-%d")
                if release_date < MIN_DATE and \
                   parts[11] in {"Complete Genome", "Chromosome", "Scaffold"} and \
                   parts[10] == "latest":
                    ftp_paths.append(parts[ftp_col_index])
            except ValueError:
                continue

print(f"[INFO] Filtered to {len(ftp_paths)} genomes released after {MIN_DATE.date()}.")

# Save ftp_paths
ftp_paths_file = METADATA_DIR / f"{CATEGORY}_ftp_paths.txt"
with open(ftp_paths_file, "w") as f:
    for path in ftp_paths:
        f.write(path + "\n")

# Construct genomic URLs
genomic_urls_file = METADATA_DIR / f"{CATEGORY}_genomic_urls.txt"
with open(genomic_urls_file, "w") as out_f:
    for path in ftp_paths:
        accession = path.split("/")[-1]
        url = f"{path}/{accession}_genomic.fna.gz"
        out_f.write(url + "\n")

# Helper function to download a file
def download_file(url, dest_path):
    try:
        with requests.get(url, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(dest_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        return True
    except Exception:
        return False

def download_worker(url):
    url = url.strip()
    filename = url.split("/")[-1]
    dest_path = GENOME_DIR / filename
    if dest_path.exists():
        return 'skipped'
    if download_file(url.replace("ftp://", "https://"), dest_path):
        return 'downloaded'
    return 'failed'

# Download genomes in parallel
print(f"[INFO] Downloading genomes for '{CATEGORY}' with {MAX_WORKERS} threads...")
downloaded = skipped = failed = 0
with open(genomic_urls_file) as f:
    urls = [line.strip() for line in f]

with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
    futures = {executor.submit(download_worker, url): url for url in urls}
    for future in as_completed(futures):
        result = future.result()
        if result == 'downloaded':
            downloaded += 1
        elif result == 'skipped':
            skipped += 1
        else:
            failed += 1

print(f"[INFO] Download complete: {downloaded} downloaded, {skipped} skipped, {failed} failed.")

# Decompress all .gz files in parallel
def decompress_gz(gz_file):
    try:
        with gzip.open(gz_file, 'rb') as f_in:
            with open(gz_file.with_suffix(''), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        gz_file.unlink()
        return True
    except Exception:
        print(f"[WARNING] Failed to decompress: {gz_file.name}")
        return False

print("[INFO] Decompressing genome files in parallel...")
gz_files = list(GENOME_DIR.glob("*.gz"))
with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
    list(executor.map(decompress_gz, gz_files))

print(f"[INFO] Done for category: {CATEGORY}")