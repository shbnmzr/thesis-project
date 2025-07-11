import argparse
import logging
from pathlib import Path


def setup_logger(log_file: Path):
    log_file.parent.mkdir(parents=True, exist_ok=True)  # Ensure log directory exists
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )


def combine_fasta(split_dir: Path, output_file: Path):
    # List of categories to expect
    categories = ["archaea", "bacteria", "fungi", "protozoa", "plasmid", "viral"]

    fasta_files = []
    for category in categories:
        cleaned_dir = split_dir / category / "cleaned"
        if cleaned_dir.exists():
            found_files = list(cleaned_dir.glob("*.fna"))
            logging.info(f"Found {len(found_files)} files in {cleaned_dir}")
            fasta_files.extend(found_files)
        else:
            logging.warning(f"Directory {cleaned_dir} does not exist, skipping.")

    if not fasta_files:
        logging.error(f"No .fna files found under {split_dir}")
        raise FileNotFoundError(f"No .fna files found under {split_dir}")

    with open(output_file, "w") as outfile:
        for fasta in sorted(fasta_files):
            logging.info(f"Adding {fasta}")
            with open(fasta, "r") as infile:
                outfile.write(infile.read())

    logging.info(f"Combined {len(fasta_files)} FASTA files into {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Combine .fna files into a single FASTA file.")
    parser.add_argument("--dataset_dir", required=True, help="Root of the dataset directory (should point to the split directory).")
    parser.add_argument("--output_file", required=True, help="Path for the combined FASTA file.")
    parser.add_argument("--log_file", default="combine.log", help="Path to the log file.")
    args = parser.parse_args()

    setup_logger(Path(args.log_file))

    split_dir = Path(args.dataset_dir).resolve()
    output_file = Path(args.output_file).resolve()

    logging.info(f"Combining genomes from {split_dir} into {output_file}")
    combine_fasta(split_dir, output_file)


if __name__ == "__main__":
    main()
