#!/bin/bash

# Step 1: Collect all FASTA files
find ./data/raw/genomes -path "./data/raw/genomes/all_genomes" -prune -o -name "*.fna" > ./data/raw/genomes/all_genomic_files.txt

# Step 2: Combine all FASTA files into one
mkdir -p ./data/raw/genomes/all_genomes
find ./data/raw/genomes -path "./data/raw/genomes/all_genomes" -prune -o -name "*.fna" -exec cp {} ./data/raw/genomes/all_genomes/ \;
cat ./data/raw/genomes/all_genomes/*.fna > ./data/raw/genomes/combined_genomes.fasta

# Check if combined FASTA is empty
if [ ! -s ./data/raw/genomes/combined_genomes.fasta ]; then
    echo "Error: Combined FASTA file is empty."
    exit 1
fi

# Step 3: Create a BLAST database
makeblastdb -in ./data/raw/genomes/combined_genomes.fasta -dbtype nucl -out custom_db

# Step 4: Perform Self-BLAST
blastn -query ./data/raw/genomes/combined_genomes.fasta -db custom_db -out ./data/raw/genomes/self_blast_results.txt \
-outfmt "6 qseqid sseqid pident length evalue bitscore" \
-perc_identity 90 -qcov_hsp_perc 80

# Check if BLAST results exist
if [ ! -f "./data/raw/genomes/self_blast_results.txt" ]; then
    echo "Error: BLAST results file './data/raw/genomes/self_blast_results.txt' not found."
    exit 1
fi

# Step 5: Filter self-matches with high similarity
awk '$1 != $2 && $3 >= 90 && $4 >= 1000 {print $1}' ./data/raw/genomes/self_blast_results.txt > ./data/raw/genomes/similar_ids.txt

# Check if seqkit is installed
if ! command -v seqkit &> /dev/null; then
    echo "Error: SeqKit is not installed. Please install it and try again."
    exit 1
fi

# Step 6: Remove highly similar sequences
seqkit grep -v -f ./data/raw/genomes/similar_ids.txt ./data/raw/genomes/combined_genomes.fasta > ./data/raw/genomes/unique_sequences.fasta

echo "Filtering complete. Unique sequences saved to 'unique_sequences.fasta'."
