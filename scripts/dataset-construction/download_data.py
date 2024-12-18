import os
from subprocess import call
import gzip
import shutil


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
            "--parallel", "10",  # Use multiple threads for faster download
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

    # Loop through all files in the input directory
    for root, dirs, files in os.walk(input_directory):
        for filename in files:
            if filename.endswith('.fna.gz'):
                # Construct the full file path
                filepath = os.path.join(root, filename)
                # Output file path by replacing .gz with .fna
                output_file = os.path.join(output_directory, filename[:-3])

                # Decompress the file
                with gzip.open(filepath, 'rb') as f_in:
                    with open(output_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

                print(f"Decompressed {filename} to {output_file}")


if __name__ == "__main__":
    # Specify the output directory for the dataset
    output_directory = "../../data/raw"  # Path where the genomes will be saved
    decompressed_output_directory = "../../data/decompressed"  # Path to store decompressed .fna files

    # Specify organism groups to download
    organism_groups_to_download = ["bacteria", "viral", "fungi", "protozoa"]

    # Step 1: Download genomes
    download_ncbi_genomes(output_directory, organism_groups_to_download)

    # Step 2: Decompress the downloaded .fna.gz files
    decompress_fna_files(output_directory, decompressed_output_directory)

    print("All downloads and decompressions completed!")
