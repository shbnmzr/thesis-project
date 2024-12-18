import os
from collections import defaultdict


def parse_md5sums(md5_file):
    """
    Parse an MD5SUMS file to extract hashes and file paths for .fna files only.

    Parameters:
    - md5_file (str): Path to the MD5SUMS file.

    Returns:
    - dict: A dictionary mapping .fna file paths to their MD5 hashes.
    """
    md5_dict = {}
    with open(md5_file, 'r') as f:
        for line in f:
            hash_value, file_path = line.strip().split("  ", 1)
            # Filter only .fna.gz files
            if file_path.endswith('.fna.gz'):
                md5_dict[file_path] = hash_value
    return md5_dict


def find_duplicate_fna_files(md5sums_dir):
    """
    Identify duplicate .fna files using MD5SUMS files in a directory.

    Parameters:
    - md5sums_dir (str): Path to the directory containing MD5SUMS files.

    Returns:
    - dict: A dictionary mapping MD5 hashes to lists of .fna file paths.
    """
    hash_to_files = defaultdict(list)

    # Traverse all subdirectories for MD5SUMS files
    for root, _, files in os.walk(md5sums_dir):
        for file in files:
            if file == "MD5SUMS":
                md5_file_path = os.path.join(root, file)
                md5_data = parse_md5sums(md5_file_path)

                # Map each hash to its associated .fna file paths
                for relative_path, hash_value in md5_data.items():
                    full_path = os.path.join(root, relative_path.strip('./'))
                    hash_to_files[hash_value].append(full_path)

    # Filter out hashes with only one file (not duplicates)
    duplicates = {hash_val: paths for hash_val, paths in hash_to_files.items() if len(paths) > 1}
    return duplicates


def main():
    # Directory containing the genomes with MD5SUMS files
    md5sums_directory = "../../data/raw/refseq"

    print(f"Scanning directory: {md5sums_directory} for duplicate .fna files using MD5SUMS...")
    duplicates = find_duplicate_fna_files(md5sums_directory)

    # Print the results
    if duplicates:
        print("\nDuplicate .fna files found:")
        for hash_value, files in duplicates.items():
            print(f"MD5: {hash_value}")
            for file in files:
                print(f"  - {file}")
            print()

        # Optional: Confirm and delete duplicates
        confirm = input("Do you want to delete duplicates, keeping only one copy? (y/n): ").strip().lower()
        if confirm == 'y':
            for hash_value, files in duplicates.items():
                for duplicate_file in files[1:]:
                    os.remove(duplicate_file)
                    print(f"Deleted duplicate: {duplicate_file}")
            print("\nAll duplicates deleted.")
    else:
        print("\nNo duplicate .fna files found.")


if __name__ == "__main__":
    main()
