import os
import gzip
import shutil
from subprocess import call
import urllib.request

import requests
from tqdm import tqdm


def download_ncbi_genomes(output_dir, organism_groups):
    """
    Downloads complete genomes of specified organism groups using ncbi-genome-download.

    Parameters:
    - output_dir (str): Directory to save the downloaded genomes.
    - organism_groups (list of str): Organism groups to download (e.g., 'bacteria', 'viral').

    Requirements:
    - ncbi-genome-download Python package (install via pip).
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for group in organism_groups:
        print(f"Downloading {group} genomes...")
        call([
            "ncbi-genome-download",
            group,
            "--formats", "fasta",
            "--assembly-level", "complete",
            "--parallel", "10",
            "--output-folder", output_dir
        ])
        print(f"Finished downloading {group} genomes!")


def decompress_fna_files(input_directory, output_directory):
    """
    Decompresses all .fna.gz files in the input directory and saves them in the output directory.

    Parameters:
    - input_directory (str): Path to the directory containing the .fna.gz files.
    - output_directory (str): Path to the directory where the decompressed .fna files will be saved.
    """
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Loop through all directories in the input directory
    for root, dirs, files in os.walk(input_directory):
        for file in files:
            # Check for .fna.gz files
            if file.endswith('.fna.gz'):
                # Full path to the .fna.gz file
                filepath = os.path.join(root, file)
                # Output file path (removing the .gz extension)
                output_file = os.path.join(output_directory, file[:-3])  # Remove .gz
                # Decompress the .fna.gz file
                with gzip.open(filepath, 'rb') as f_in:
                    with open(output_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)


def extract_fna_files_from_md5sum(input_directory, output_directory):
    """
    Processes MD5SUMS files to extract FTP links for .fna.gz files and downloads them, avoiding redundant downloads.

    Parameters:
    - input_directory (str): Directory containing subdirectories with MD5SUMS files.
    - output_directory (str): Directory to store the downloaded .fna files.
    """

    def download_file_with_progress(url, output_file):
        if os.path.exists(output_file):
            print(f"File {output_file} already exists, skipping download.")
            return

        response = requests.get(url, stream=True)
        total_size_in_bytes= int(response.headers.get('content-length', 0))
        block_size = 1024 # 1 Kibibyte
        with open(output_file, 'wb') as file:
            for data in tqdm(response.iter_content(block_size), total=total_size_in_bytes/block_size, unit='KB', unit_scale=True, desc=url):
                file.write(data)

    # Loop through each subdirectory in the input directory
    for root, dirs, files in os.walk(input_directory):
        if 'MD5SUMS' in files:
            md5sum_file_path = os.path.join(root, 'MD5SUMS')
            print(f"Processing MD5SUMS in {root}...")

            # Read the MD5SUMS file to get the FTP links for .fna.gz files
            with open(md5sum_file_path, 'r') as md5_file:
                for line in tqdm(md5_file, desc="Processing MD5SUMS"):
                    # Extract FTP link for .fna.gz files
                    parts = line.split()
                    if len(parts) >= 2 and parts[1].endswith('.fna.gz'):
                        # Full FTP link to the .fna.gz file
                        ftp_link = parts[1][2:]  # Remove './' prefix from the link

                        # Construct the full output file path
                        output_file = os.path.join(output_directory, parts[1][2:])

                        # Download the file using requests with progress bar, avoiding redundant downloads
                        download_file_with_progress(ftp_link, output_file)


def main():
    # Specify the output directory for the dataset
    output_directory = "../../data/raw"

    # Specify organism groups to download
    # Groups: 'bacteria', 'archaea', 'viral', 'fungi', 'protozoa', etc.
    organism_groups_to_download = ["bacteria", "viral", "fungi", "protozoa"]

    # Start the download process
    download_ncbi_genomes(output_directory, organism_groups_to_download)

    # Specify the directory containing the downloaded data
    input_directory = "../../data/raw/refseq"
    # Specify the output directory for decompressed .fna files
    decompressed_output_directory = "../../data/decompressed"

    decompress_fna_files(input_directory, decompressed_output_directory)
    # Extract and decompress .fna.gz files
    # extract_fna_files_from_md5sum(input_directory, decompressed_output_directory)

    print("All downloads and decompressions completed!")


if __name__ == "__main__":
    main()
