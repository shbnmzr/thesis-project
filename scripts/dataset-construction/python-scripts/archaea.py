import os
import requests
import time
import gzip
import shutil
from pathlib import Path

CATEGORY = "archaea"
BASE_DIR = Path("./4CAC_dataset")
METADATA_DIR = BASE_DIR / "metadata"
GENOME_DIR = BASE_DIR / "genomes" / CATEGORY

# Create necessary directories
GENOME_DIR.mkdir(parents=True, exist_ok=True)
METADATA_DIR.mkdir(parents=True, exist_ok=True)

# Step 1: Download the assembly summary
assembly_summary_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/genbank/{CATEGORY}/assembly_summary.txt"
summary_file = METADATA_DIR / f"{CATEGORY}_assembly_summary.txt"

print(f"Downloading assembly summary for {CATEGORY}...", flush=True)
response = requests.get(assembly_summary_url)
response.raise_for_status()
with open(summary_file, "w") as f:
    f.write(response.text)

# Find ftp_path column index
with open(summary_file) as f:
    header = next(line for line in f if line.startswith("#assembly_accession")).strip().split('\t')
ftp_col_index = header.index("ftp_path")

# Extract relevant ftp paths
ftp_paths = []
with open(summary_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split('\t')
        if len(parts) > 11 and parts[11] in {"Complete Genome", "Chromosome", "Scaffold"} and parts[10] == "latest":
            ftp_paths.append(parts[ftp_col_index])

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
    except Exception as e:
        print(f"Failed to download {url}: {e}", flush=True)
        return False

# Download genomes
print(f"Downloading genomes for {CATEGORY}...", flush=True)
with open(genomic_urls_file) as f:
    for url in f:
        url = url.strip()
        filename = url.split("/")[-1]
        dest_path = GENOME_DIR / filename
        if not dest_path.exists():
            print(f"Downloading {filename}...", flush=True)
            success = download_file(url.replace("ftp://", "https://"), dest_path)
            if success:
                print(f"{filename} downloaded successfully.", flush=True)
            else:
                print(f"Skipping {filename} due to error.", flush=True)
            time.sleep(1)
        else:
            print(f"{filename} already downloaded. Skipping.", flush=True)

# Decompress all .gz files
for gz_file in GENOME_DIR.glob("*.gz"):
    try:
        with gzip.open(gz_file, 'rb') as f_in:
            with open(gz_file.with_suffix(''), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        gz_file.unlink()  # Remove the .gz file after decompression
    except Exception as e:
        print(f"Failed to decompress {gz_file.name}: {e}", flush=True)

print(f"Completed downloading for {CATEGORY}.", flush=True)