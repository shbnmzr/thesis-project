import argparse
import subprocess
import logging
from pathlib import Path


def setup_logger(log_file: Path):
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )


def simulate_reads(genomes_fasta: Path, output_prefix: Path, n_reads: int, cpus: int):
    cmd = [
        "iss", "generate",
        "--genomes", str(genomes_fasta),
        "--model", "miseq",
        "--output", str(output_prefix),
        "--n_reads", str(n_reads),
        "--cpus", str(cpus)
    ]

    logging.info(f"Running InSilicoSeq with command: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
        logging.info(f"InSilicoSeq completed successfully. Reads generated at {output_prefix}_R1.fastq and {output_prefix}_R2.fastq")
    except subprocess.CalledProcessError as e:
        logging.error(f"InSilicoSeq failed with error: {e}")
        raise


def main():
    parser = argparse.ArgumentParser(description="Simulate Illumina MiSeq reads using InSilicoSeq.")
    parser.add_argument("--genomes_fasta", required=True, help="Path to the combined genomes FASTA file.")
    parser.add_argument("--output_dir", required=True, help="Directory to save the simulated reads.")
    parser.add_argument("--output_prefix", default="simulated", help="Prefix for the output FASTQ files.")
    parser.add_argument("--n_reads", type=int, default=10000000, help="Total number of reads to simulate.")
    parser.add_argument("--cpus", type=int, default=8, help="Number of CPU threads to use.")
    parser.add_argument("--log_file", default="simulate_reads.log", help="Path to the log file.")
    args = parser.parse_args()

    setup_logger(Path(args.log_file))

    genomes_fasta = Path(args.genomes_fasta).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    output_prefix = output_dir / args.output_prefix

    logging.info(f"Starting read simulation from {genomes_fasta} into {output_dir} with prefix {args.output_prefix}")

    simulate_reads(genomes_fasta, output_prefix, args.n_reads, args.cpus)


if __name__ == "__main__":
    main()
