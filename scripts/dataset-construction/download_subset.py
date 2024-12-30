import os
import pandas as pd
import subprocess
import zipfile

# Number of records to download per metadata file (for local testing).
NUM_RECORDS = 5


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

        # ------------------------------------------------------------
        # For local testing, limit the number of records downloaded:
        df = df.head(NUM_RECORDS)
        # ------------------------------------------------------------

        # Output folder for the taxonomic group
        group_folder = os.path.join(output_dir, group)

        # Download each genome
        for accession in df["Accession"]:
            # Skip download if there's already at least one .fna file for this accession
            if already_downloaded(accession, group_folder):
                print(f"Skipping {accession}: already downloaded.")
                continue

            download_genome(accession, group_folder)

    print("Download process completed.")


# Function to check if the .fna files for a given accession might already exist
def already_downloaded(accession, group_folder):
    """
    Check if at least one .fna file with the accession in its name
    already exists in group_folder.
    """
    for fname in os.listdir(group_folder):
        if fname.startswith(accession) and fname.endswith(".fna"):
            return True
    return False


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

        # Delete the ZIP file after extraction
        os.remove(zip_file)
        print(f"Deleted ZIP file: {zip_file}")
    except subprocess.CalledProcessError as e:
        print(f"Failed to download {accession}: {e}")
    except Exception as e:
        print(f"An error occurred while processing {accession}: {e}")


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
