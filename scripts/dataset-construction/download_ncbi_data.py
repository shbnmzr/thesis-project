import os
import pandas as pd
import subprocess
import zipfile


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

    # Process each metadata file
    for group, path in file_paths.items():
        print(f"Processing {group}...")

        # Load the metadata file
        df = pd.read_csv(path, sep="\t", header=None, names=["Accession", "Release Date"])

        # Output folder for the taxonomic group
        group_folder = os.path.join(output_dir, group)

        # Download each genome
        for accession in df["Accession"]:
            download_genome(accession, group_folder)

    print("Download process completed.")


# Function to download genomes
def download_genome(accession, output_folder):
    try:
        # Define the output file path
        zip_file = os.path.join(output_folder, f"{accession}.zip")

        # Command to download genome data using `datasets`
        command = [
            "datasets",
            "download",
            "genome",
            "accession",
            accession,
            "--filename",
            zip_file,
        ]

        # Execute the command
        subprocess.run(command, check=True)
        print(f"Downloaded: {accession} to {zip_file}")

        # Extract .fna files from the ZIP archive
        extract_fna(zip_file, output_folder)
    except subprocess.CalledProcessError as e:
        print(f"Failed to download {accession}: {e}")


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
                print(f"Extracted .fna files: {fna_files} to {output_folder}")
            else:
                print(f"No .fna files found in {zip_file}")
    except zipfile.BadZipFile as e:
        print(f"Failed to extract {zip_file}: {e}")


if __name__ == "__main__":
    main()
