#!/bin/bash
set -e

CATEGORY="protozoa"
BASE_DIR="/home/zareshahra/thesis/4CAC_dataset"
METADATA_DIR="$BASE_DIR/metadata"
GENOME_DIR="$BASE_DIR/genomes/$CATEGORY"

mkdir -p "$GENOME_DIR" "$METADATA_DIR"

# STEP 1: Clearly download the assembly summary
curl -s "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/$CATEGORY/assembly_summary.txt" \
    -o "$METADATA_DIR/${CATEGORY}_assembly_summary.txt"

summary_file="$METADATA_DIR/${CATEGORY}_assembly_summary.txt"

ftp_col=$(grep "^#assembly_accession" "$summary_file" | head -n 1 | \
    awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="ftp_path") print i; exit}')

awk -F '\t' -v ftp_col="$ftp_col" \
    'NR>1 && ($12=="Complete Genome" || $12=="Chromosome" || $12=="Scaffold") && $11=="latest" {print $ftp_col}' \
    "$summary_file" > "$METADATA_DIR/${CATEGORY}_ftp_paths.txt"

awk '{print $1 "/" gensub(".*/","","g",$1)"_genomic.fna.gz"}' \
    "$METADATA_DIR/${CATEGORY}_ftp_paths.txt" > "$METADATA_DIR/${CATEGORY}_genomic_urls.txt"

echo "Downloading genomes for $CATEGORY..."
while read -r url; do
    filename=$(basename "$url")
    if [ ! -f "$GENOME_DIR/$filename" ]; then
        wget -q -nc "$url" -P "$GENOME_DIR"
    else
        echo "$filename already downloaded. Skipping."
    fi
done < "$METADATA_DIR/${CATEGORY}_genomic_urls.txt"

gunzip -f "$GENOME_DIR"/*.gz 2>/dev/null || true

echo "Completed downloading for $CATEGORY."
