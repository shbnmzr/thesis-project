import gzip
import os
import shutil
from subprocess import call


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


def count_compressed_files(input_directory):
    """
    Counts the number of .fna.gz files in the input directory for each organism group.

    Parameters:
    - input_directory (str): Path to the directory containing .fna.gz files.

    Returns:
    - dict: A dictionary where keys are organism groups and values are counts of .fna.gz files.
    """
    counts = {}
    for root, _, files in os.walk(input_directory):
        # Use the subdirectory name to determine the organism group
        group = os.path.basename(root)
        fna_gz_count = sum(1 for file in files if file.endswith('.fna.gz'))
        if fna_gz_count > 0:
            counts[group] = counts.get(group, 0) + fna_gz_count
    return counts


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


def main():
    # Specify the output directory for the dataset
    output_directory = "../../data/raw"

    # Specify organism groups to download
    organism_groups_to_download = ["bacteria", "viral", "fungi", "protozoa"]

    # Start the download process
    download_ncbi_genomes(output_directory, organism_groups_to_download)

    # Specify the directory containing the downloaded data
    input_directory = "../../data/raw/refseq"

    # Count the number of compressed files (.fna.gz) for each group
    compressed_file_counts = count_compressed_files(input_directory)

    # Print the counts
    print("\nCounts of compressed .fna.gz files per organism group:")
    for group, count in compressed_file_counts.items():
        print(f"{group}: {count} files")

    print("\nAll downloads and counting completed!")


if __name__ == "__main__":
    main()
