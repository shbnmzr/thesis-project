#!/usr/bin/env python3
"""
Download and prepare genomes from NCBI GenBank for a given category.

Usage:
    python download_genomes.py --category archaea
    python download_genomes.py --category bacteria --base_dir ./my_dataset
"""

import argparse
import requests
import time
import gzip
import shutil
import logging
from pathlib import Path
from typing import List

# -----------------------------
# CONFIGURATION
# -----------------------------
MAX_RETRIES = 3
RETRY_DELAY = 5

# -----------------------------
# LOGGING CONFIGURATION
# -----------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler()]
)


# -----------------------------
# UTILITY FUNCTIONS
# -----------------------------
def ensure_directories(genome_dir: Path, metadata_dir: Path) -> None:
    """
    Create necessary directories.
    """
    genome_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)


def download_text_file(url: str, dest: Path) -> None:
    """
    Download a text file and save it locally.
    """
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            logging.info(f"Downloading: {url}")
            response = requests.get(url, timeout=60)
            response.raise_for_status()
            dest.write_text(response.text)
            return
        except Exception as e:
            logging.warning(f"Attempt {attempt} failed for {url}: {e}")
            time.sleep(RETRY_DELAY)
    raise RuntimeError(f"Failed to download {url} after {MAX_RETRIES} attempts.")


def parse_ftp_paths(summary_file: Path) -> List[str]:
    """
    Extract valid FTP paths from assembly summary.
    """
    with summary_file.open() as f:
        header = next(line for line in f if line.startswith("#assembly_accession")).strip().split("\t")
    ftp_col_index = header.index("ftp_path")

    ftp_paths = []
    with summary_file.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) > 11 and parts[11] in {"Complete Genome", "Chromosome", "Scaffold"} and parts[10] == "latest":
                ftp_paths.append(parts[ftp_col_index])
    return ftp_paths


def save_list_to_file(data: List[str], dest: Path) -> None:
    """
    Save a list of strings to a file.
    """
    dest.write_text("\n".join(data))


def download_file(url: str, dest_path: Path) -> bool:
    """
    Download a file with retries.
    """
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            with requests.get(url, stream=True, timeout=60) as r:
                r.raise_for_status()
                with open(dest_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
            return True
        except Exception as e:
            logging.warning(f"Attempt {attempt} failed for {url}: {e}")
            time.sleep(RETRY_DELAY)
    logging.error(f"Failed to download {url} after {MAX_RETRIES} attempts.")
    return False


def decompress_gz_files(directory: Path) -> None:
    """
    Decompress all .gz files in a directory.
    """
    for gz_file in directory.glob("*.gz"):
        try:
            with gzip.open(gz_file, "rb") as f_in, open(gz_file.with_suffix(""), "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            gz_file.unlink()
            logging.info(f"Decompressed {gz_file.name}")
        except Exception as e:
            logging.error(f"Failed to decompress {gz_file.name}: {e}")


def main(category: str, base_dir: Path):
    """
    Main pipeline execution.
    """
    metadata_dir = base_dir / "metadata"
    genome_dir = base_dir / "genomes" / category
    ensure_directories(genome_dir, metadata_dir)

    # Step 1: Download metadata
    summary_file = metadata_dir / f"{category}_assembly_summary.txt"
    url = f"https://ftp.ncbi.nlm.nih.gov/genomes/genbank/{category}/assembly_summary.txt"
    download_text_file(url, summary_file)

    # Step 2: Extract and save FTP paths
    ftp_paths = parse_ftp_paths(summary_file)
    save_list_to_file(ftp_paths, metadata_dir / f"{category}_ftp_paths.txt")

    # Step 3: Construct genomic URLs
    genomic_urls = [f"{p}/{p.split('/')[-1]}_genomic.fna.gz" for p in ftp_paths]
    save_list_to_file(genomic_urls, metadata_dir / f"{category}_genomic_urls.txt")

    # Step 4: Download genomes
    for url in genomic_urls:
        filename = url.split("/")[-1]
        dest_path = genome_dir / filename
        if dest_path.exists():
            logging.info(f"{filename} already exists. Skipping.")
            continue
        https_url = url.replace("ftp://", "https://")
        if download_file(https_url, dest_path):
            logging.info(f"{filename} downloaded successfully.")
        else:
            logging.error(f"Skipping {filename} due to repeated failures.")
        time.sleep(1)

    # Step 5: Decompress genomes
    decompress_gz_files(genome_dir)
    logging.info(f"Completed processing for {category}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download genomes from NCBI GenBank.")
    parser.add_argument("--category", required=True, help="Category (e.g., archaea, bacteria, viral)")
    parser.add_argument("--base_dir", default="./4CAC_dataset", help="Base directory for storing data")
    args = parser.parse_args()

    main(args.category, Path(args.base_dir))
