#!/bin/bash
set -e

BASE_DIR="./dataset"
CATEGORY="plasmids"
GENOME_DIR="$BASE_DIR/genomes/$CATEGORY"
METADATA_DIR="$BASE_DIR/metadata"

mkdir -p "$GENOME_DIR" "$METADATA_DIR"

echo "Downloading all plasmid genomes (RefSeq)..."

curl -s "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/" | \
grep -o "plasmid.*genomic.fna.gz" | \
sed 's|^|ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/|' \
> "$METADATA_DIR/plasmid_file_urls.txt"

while read -r url; do
    filename=$(basename "$url")
    if [ ! -f "$GENOME_DIR/$filename" ]; then
        wget -q -nc "$url" -P "$GENOME_DIR"
    else
        echo "$filename already downloaded. Skipping."
    fi
done < "$METADATA_DIR/plasmid_file_urls.txt"

gunzip -f "$GENOME_DIR"/*.gz 2>/dev/null || true
cat "$GENOME_DIR"/*.fna > "$GENOME_DIR/all_plasmids.fasta"

echo "Plasmid genomes downloaded successfully."
