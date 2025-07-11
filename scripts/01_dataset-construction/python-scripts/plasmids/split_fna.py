from pathlib import Path
from Bio import SeqIO
import os

INPUT_FASTA = Path("../4CAC_dataset/genomes/plasmids/all_plasmids.fasta")
OUTPUT_DIR = Path("../4CAC_dataset/train/plasmid/raw")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

genomes = {}

for record in SeqIO.parse(INPUT_FASTA, "fasta"):
    header = record.id
    accession = header.split("|")[0]  # e.g., GCF_XXXXX
    if accession not in genomes:
        genomes[accession] = []
    genomes[accession].append(record)

print(f"[INFO] Found {len(genomes)} unique plasmid accessions.")

for accession, records in genomes.items():
    output_path = OUTPUT_DIR / f"{accession}.fna"
    SeqIO.write(records, output_path, "fasta")

print(f"[DONE] Saved split FASTA files to: {OUTPUT_DIR}")