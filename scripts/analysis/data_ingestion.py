# scripts/data_ingestion.py
from Bio import SeqIO
import os


def stream_fasta_files(directory):
    """
    Recursively yields each FASTA filename + the SeqIO parse iterator.
    This doesn't store sequences in memory; it just creates a generator
    for each file ending with .fna, .fa, or .fasta.
    """
    for root, dirs, files in os.walk(directory):
        for fname in files:
            if fname.lower().endswith(('.fna', '.fa', '.fasta')):
                file_path = os.path.join(root, fname)
                yield file_path, SeqIO.parse(file_path, "fasta")


def main():
    data_dir = '../../data/raw/genomes'

    file_count = 0
    seq_count = 0

    for file_path, seq_generator in stream_fasta_files(data_dir):
        file_count += 1
        # Process or count sequences on the fly
        for rec in seq_generator:
            seq_count += 1

    print("Total files:", file_count)
    print("Total sequences:", seq_count)


if __name__ == '__main__':
    main()
