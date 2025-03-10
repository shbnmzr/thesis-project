#!/bin/bash
set -e

CATEGORY="bacteria"
BASE_DIR="/home/zareshahra/thesis/4CAC_dataset"
METADATA_DIR="$BASE_DIR/metadata"
GENOME_DIR="$BASE_DIR/genomes/$CATEGORY"

mkdir -p "$GENOME_DIR" "$METADATA_DIR"

# Step 1: Download bacteria assembly summary from GenBank
if [ ! -f "$METADATA_DIR/${CATEGORY}_assembly_summary.txt" ]; then
    echo "Downloading assembly summary for $CATEGORY..."
    curl -s "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/$CATEGORY/assembly_summary.txt" \
         -o "$METADATA_DIR/${CATEGORY}_assembly_summary.txt"
else
    echo "Assembly summary already downloaded. Skipping."
fi

summary_file="$METADATA_DIR/${CATEGORY}_assembly_summary.txt"

# Step 2: Dynamically find the FTP path column
ftp_col=$(grep "^#assembly_accession" "$summary_file" | head -n 1 | \
    awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="ftp_path") print i; exit}')

# Verify ftp_col found
if [ -z "$ftp_col" ]; then
    echo "FTP path column not found! Exiting." >&2
    exit 1
fi

# Step 3: Extract FTP paths for complete genomes, latest only
awk -F '\t' -v ftp_col="$ftp_col" \
    'NR>1 && $12=="Complete Genome" && $11=="latest" && $ftp_col!="na" {print $ftp_col}' \
    "$summary_file" > "$METADATA_DIR/${CATEGORY}_ftp_paths.txt"

# Step 4: Convert to direct genomic FASTA URLs
awk '{print $1 "/" gensub(".*/","","g",$1)"_genomic.fna.gz"}' \
    "$METADATA_DIR/${CATEGORY}_ftp_paths.txt" > "$METADATA_DIR/${CATEGORY}_genomic_urls.txt"

# Step 4: Download genomes explicitly checking for existing files
while read -r url; do
    filename=$(basename "$url")
    if [ ! -f "$GENOME_DIR/$filename" ]; then
        wget -q -nc "$url" -P "$GENOME_DIR"
    else
        echo "$filename already downloaded. Skipping."
    fi
done < "$METADATA_DIR/${CATEGORY}_genomic_urls.txt"

# Step 5: Extract genome files
gunzip -f "$GENOME_DIR"/*.gz 2>/dev/null || true

# Step 6: Verify the total downloaded sequences
cat "$GENOME_DIR"/*.fna > "$GENOME_DIR/all_bacteria.fasta"

sequence_count=$(grep -c "^>" "$GENOME_DIR/all_bacteria.fasta" || echo 0)
echo "Completed downloading for $CATEGORY. Total sequences downloaded: $sequence_count"