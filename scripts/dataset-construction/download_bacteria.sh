# Step 1: Download the assembly summary file for the targeted group (e.g., bacteria)
curl 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt' > ./data/raw/metadata/bacterial_summary.txt

# Step 2: Filter for complete genomes (or other criteria relevant to 4CAC)
awk '{FS="\t"} !/^#/ {if($12=="Complete Genome" && $11=="latest") print $20}' ./data/raw/metadata/bacterial_summary.txt > ./data/raw/metadata/bacterial_genbank_urls.txt

# Step 3: Modify URLs to point to specific file types needed for 4CAC (e.g., genomic FASTA)
sed -r 's|(https://ftp.ncbi.nlm.nih.gov/genomes/all/.+/)(GCA_.+)|\1\2/\2_genomic.fna.gz|' ./data/raw/metadata/bacterial_genbank_urls.txt > ./data/raw/metadata/bacterial_genomic_file

# Step 4: Create a directory to store downloaded genomes
mkdir -p ./data/raw/genomes/bacterial

# Step 5: Download all genome files listed in the URL file
wget --input-file=./data/raw/metadata/bacterial_genomic_file -P ./data/raw/genomes/bacterial

# Step 6: Decompress all downloaded files
gunzip ./data/raw/genomes/bacterial/*.gz 2>/dev/null || echo "No files to decompress"

# Step 7: Combine all decompressed FASTA files into a single file
if ls ./data/raw/genomes/bacterial/*.fna 1> /dev/null 2>&1; then
    cat ./data/raw/genomes/bacterial/*.fna > ./data/raw/genomes/bacterial/all_bacterial.fasta
else
    echo "No .fna files to combine"
fi
