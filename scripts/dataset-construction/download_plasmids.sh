#!/bin/bash
set -e

BASE_DIR="./4CAC_dataset"
CATEGORY="plasmids"
GENOME_DIR="$BASE_DIR/genomes/$CATEGORY"
METADATA_DIR="$BASE_DIR/metadata"

mkdir -p "$GENOME_DIR" "$METADATA_DIR"

echo "Downloading all plasmid genomes (RefSeq)..."

# Step 1: Get URLs for genomic FASTA files explicitly from FTP directory
curl -s "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/" | \
grep -o "plasmid.*genomic.fna.gz" | \
sed 's|^|ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/|' \
> "$METADATA_DIR/plasmid_file_urls.txt"

# Step 2: Download files explicitly
while read -r url; do
    filename=$(basename "$url")
    if [ ! -f "$GENOME_DIR/$filename" ]; then
        wget -q -nc "$url" -P "$GENOME_DIR"
    else
        echo "$filename already downloaded. Skipping."
    fi
done < "$METADATA_DIR/plasmid_file_urls.txt"

# Step 3: Explicitly decompress files safely
echo "Decompressing plasmid genome files..."
find "$GENOME_DIR" -type f -name "*.gz" -exec gunzip -f {} +

# Step 4: Combine explicitly into one FASTA file
echo "Combining decompressed plasmid genomes..."
cat "$GENOME_DIR"/*.fna > "$GENOME_DIR/all_plasmids.fasta"

# Explicitly verify the combined FASTA file
sequence_count=$(grep -c "^>" "$GENOME_DIR/all_plasmids.fasta" || echo 0)
echo "Explicitly downloaded plasmid sequences: $sequence_count"

echo "Plasmid genomes downloaded, decompressed, and combined successfully!"
