#!/bin/bash
set -e

CATEGORY="viruses"
BASE_DIR="/home/zareshahra/thesis/4CAC_dataset"
METADATA_DIR="$BASE_DIR/metadata"
GENOME_DIR="$BASE_DIR/genomes/$CATEGORY"

mkdir -p "$GENOME_DIR" "$METADATA_DIR"

# Step 1: Download viral assembly summary from GenBank
if [ ! -f "$METADATA_DIR/${CATEGORY}_assembly_summary.txt" ]; then
    echo "Downloading assembly summary for $CATEGORY..."
    curl -s "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary.txt" \
         -o "$METADATA_DIR/${CATEGORY}_assembly_summary.txt"
else
    echo "Assembly summary already downloaded. Skipping."
fi

summary_file="$METADATA_DIR/${CATEGORY}_assembly_summary.txt"

# Step 2: Dynamically find the FTP path column
ftp_col=$(grep "^#assembly_accession" "$summary_file" | head -n 1 | \
    awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="ftp_path") print i; exit}')

if [ -z "$ftp_col" ]; then
    echo "Failed to find FTP path column. Exiting."
    exit 1
fi

# Step 3: Extract FTP paths for complete genomes, latest only
awk -F '\t' -v ftp_col="$ftp_col" \
    'NR>1 && ($12=="Complete Genome" || $12=="Chromosome" || $12=="Scaffold") && $11=="latest" && $ftp_col!="na" {print $ftp_col}' \
    "$METADATA_DIR/${CATEGORY}_assembly_summary.txt" > "$METADATA_DIR/${CATEGORY}_ftp_paths.txt"

# Step 4: Convert to direct URLs for genome FASTA downloads
awk '{print $1 "/" gensub(".*/", "", "g", $1) "_genomic.fna.gz"}' \
    "$METADATA_DIR/${CATEGORY}_ftp_paths.txt" > "$METADATA_DIR/${CATEGORY}_genomic_urls.txt"

# Step 5: Download genomes, skipping already downloaded files
while read -r url; do
    filename=$(basename "$url")
    if [ ! -f "$BASE_DIR/genomes/$CATEGORY/$filename" ]; then
        wget -q -nc "$url" -P "$BASE_DIR/genomes/$CATEGORY"
    else
        echo "$filename already downloaded. Skipping."
    fi
done < "$METADATA_DIR/${CATEGORY}_genomic_urls.txt"

# Step 5: Extract genome sequences
gunzip -f "$GENOME_DIR"/*.gz 2>/dev/null || true

# Step 6: Combine all genome sequences into one FASTA file
cat "$GENOME_DIR"/*.fna > "$GENOME_DIR/all_${CATEGORY}.fasta"

# Verify explicitly downloaded sequences
sequence_count=$(grep -c "^>" "$GENOME_DIR/all_${CATEGORY}.fasta" || echo 0)
echo "$CATEGORY genomes downloaded from GenBank: $sequence_count sequences"


echo "Viral genomes downloaded and combined successfully."