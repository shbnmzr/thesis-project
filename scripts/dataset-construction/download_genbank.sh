#!/bin/bash
set -e  # Exit on any error

BASE_DIR="./4CAC_dataset"
mkdir -p "$BASE_DIR/metadata"
CATEGORIES=("bacteria" "archaea" "fungi" "protozoa" "viruses" "plasmids")

# Create separate directories for each genome category
for category in "${CATEGORIES[@]}"; do
    mkdir -p "$BASE_DIR/genomes/$category"
done

# 1. Download assembly summaries for Bacteria, Archaea, Fungi, and Protozoa (GenBank for completeness)
for category in bacteria archaea fungi protozoa; do
    echo "Fetching $category assembly summary..."
    curl -s "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/${category}/assembly_summary.txt" \
         -o "$BASE_DIR/metadata/${category}_assembly_summary.txt"
done

# 2. Parse each summary to get FTP paths of complete genomes (latest versions only)
for category in bacteria archaea fungi protozoa; do
    echo "Processing $category genome list..."
    awk -F $'\t' 'NR>1 && $12=="Complete Genome" && $11=="latest" {print $20}' \
        "$BASE_DIR/metadata/${category}_assembly_summary.txt" \
        > "$BASE_DIR/metadata/${category}_ftp_paths.txt"

    # Convert base FTP path to direct URL for the genomic FASTA (.fna.gz)
    sed -r 's|(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/.+/)(GC[AF]_.+)|\1\2/\2_genomic.fna.gz|' \
        "$BASE_DIR/metadata/${category}_ftp_paths.txt" \
        > "$BASE_DIR/metadata/${category}_genomic_urls.txt"

    # Download each genome into its respective folder
    echo "Downloading genomes for $category..."
    wget -q -nc -i "$BASE_DIR/metadata/${category}_genomic_urls.txt" -P "$BASE_DIR/genomes/$category"
    
    # Extract genome sequences
    gunzip -f "$BASE_DIR/genomes/$category"/*.gz  2>/dev/null
done

# 3. Download all viral genomes from GenBank (as RefSeq lacks full virus dataset)
echo "Downloading all viral genomes (GenBank)..."
mkdir -p "$BASE_DIR/genomes/viruses"
wget -q -O "$BASE_DIR/viruses_all.fna.tar.gz" "ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/all.fna.tar.gz"
tar -xzf "$BASE_DIR/viruses_all.fna.tar.gz" -C "$BASE_DIR/genomes/viruses"
# Combine all extracted viral sequences into one file
find "$BASE_DIR/genomes/viruses" -name "*.fna" -exec cat {} + > "$BASE_DIR/genomes/viruses/all_viruses.fasta"

# 4. Download all plasmid genomes from RefSeq
echo "Downloading all plasmid genomes (RefSeq)..."
mkdir -p "$BASE_DIR/genomes/plasmids"
curl -s "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid" | grep "genomic.fna.gz" \
     | sed 's|^|ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/|' \
     > "$BASE_DIR/metadata/plasmid_file_urls.txt"
wget -q -nc -i "$BASE_DIR/metadata/plasmid_file_urls.txt" -P "$BASE_DIR/genomes/plasmids"
gunzip -f "$BASE_DIR/genomes/plasmids"/*.gz  2>/dev/null
cat "$BASE_DIR/genomes/plasmids"/*.fna > "$BASE_DIR/genomes/plasmids/all_plasmids.fasta"

# 5. Verify Downloaded Files
echo "Verifying downloaded genome counts..."
for category in "${CATEGORIES[@]}"; do
    count=$(grep -c ">" "$BASE_DIR/genomes/$category/all_${category}.fasta" 2>/dev/null || echo 0)
    echo "$category: $count sequences downloaded"
done

echo "All expected genomes downloaded. Organized in: $BASE_DIR/genomes/"