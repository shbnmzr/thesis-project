import os
import time

import pandas as pd
import subprocess
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed


def main():
    # Define metadata file paths and output directories
    file_paths = {
        "Bacteria": "../../data/raw/metadata/bacteria_ids.tsv",
        "Archaea": "../../data/raw/metadata/archaea_ids.tsv",
        "Fungi": "../../data/raw/metadata/fungi_ids.tsv",
        "Eukaryote": "../../data/raw/metadata/eukaryote_ids.tsv",
        "Viruses": "../../data/raw/metadata/viruses_ids.tsv",
        "Plasmids": "../../data/raw/metadata/plasmids_ids.tsv",
    }

    output_dir = "../../data/raw/genomes"

    # Ensure output directories exist
    os.makedirs(output_dir, exist_ok=True)
    for group in file_paths.keys():
        os.makedirs(os.path.join(output_dir, group), exist_ok=True)

    # Use ThreadPoolExecutor for parallel downloads
    with ThreadPoolExecutor() as executor:
        future_to_accession = {}

        # Process each metadata file
        for group, path in file_paths.items():
            print(f"Processing {group}...")

            # Load the metadata file
            df = pd.read_csv(path, sep="\t", header=None, names=["Accession", "Release Date"])

            # Output folder for the taxonomic group
            group_folder = os.path.join(output_dir, group)

            # Submit genome downloads to the thread pool
            for accession in df["Accession"]:
                future = executor.submit(download_genome, accession, group_folder)
                future_to_accession[future] = accession

        # Monitor progress
        for future in as_completed(future_to_accession):
            accession = future_to_accession[future]
            try:
                future.result()
                print(f"Completed: {accession}")
            except Exception as e:
                print(f"Failed to process {accession}: {e}")

    print("Download process completed.")


# Function to download genomes
def download_genome(accession, output_folder, retries=3):
    for attempt in range(retries):
        try:
            fna_file = os.path.join(output_folder, f"{accession}.fna")
            if os.path.exists(fna_file):
                print(f"Skipped: {accession} (already downloaded)")
                return

            zip_file = os.path.join(output_folder, f"{accession}.zip")
            command = [
                "datasets", "download", "genome", "accession",
                accession, "--filename", zip_file
            ]
            subprocess.run(command, check=True)
            extract_fna(zip_file, output_folder)
            os.remove(zip_file)
            return
        except subprocess.CalledProcessError as e:
            print(f"Attempt {attempt + 1} failed for {accession}: {e}")
            if attempt < retries - 1:
                time.sleep(5)
            else:
                print(f"Failed after {retries} attempts: {accession}")


# Function to extract .fna files from a ZIP archive
def extract_fna(zip_file, output_folder):
    try:
        with zipfile.ZipFile(zip_file, 'r') as zip_ref:
            # List all files in the ZIP archive
            files = zip_ref.namelist()

            # Extract only .fna files
            fna_files = [f for f in files if f.endswith(".fna")]
            zip_ref.extractall(output_folder, members=fna_files)

            if fna_files:
                # Rename the extracted .fna file to match the accession
                extracted_file = os.path.join(output_folder, fna_files[0])
                renamed_file = os.path.join(output_folder, f"{os.path.basename(zip_file).replace('.zip', '.fna')}")
                os.rename(extracted_file, renamed_file)
                print(f"Extracted and renamed: {renamed_file}")
            else:
                print(f"No .fna files found in {zip_file}")
    except zipfile.BadZipFile as e:
        print(f"Failed to extract {zip_file}: {e}")


if __name__ == "__main__":
    main()
