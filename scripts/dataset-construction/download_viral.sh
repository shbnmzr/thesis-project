# Step 1: Download the assembly summary file for the targeted group
curl 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt' > ./data/raw/metadata/viral_summary.txt

# Step 2: Filter for complete genomes (or other criteria relevant to 4CAC)
awk '{FS="\t"} !/^#/ {if($12=="Complete Genome" && $11=="latest") print $20}' ./data/raw/metadata/viral_summary.txt > ./data/raw/metadata/viral_genbank_urls.txt

# Step 3: Modify URLs to point to specific file types needed for 4CAC (e.g., genomic FASTA)
sed -r 's|(https://ftp.ncbi.nlm.nih.gov/genomes/all/.+/)(GCA_.+)|\1\2/\2_genomic.fna.gz|' ./data/raw/metadata/viral_genbank_urls.txt > ./data/raw/metadata/viral_genomic_file

# Step 4: Create a directory to store downloaded genomes
mkdir -p ./data/raw/genomes/viral

# Step 5: Download all genome files listed in the URL file
wget --input-file=./data/raw/metadata/viral_genomic_file -P ./data/raw/genomes/viral

# Step 6: Decompress all downloaded files
gunzip ./data/raw/genomes/viral/*.gz 2>/dev/null || echo "No files to decompress"

# Step 7: Combine all decompressed FASTA files into a single file
if ls ./data/raw/genomes/viral/*.fna 1> /dev/null 2>&1; then
    cat ./data/raw/genomes/viral/*.fna > ./data/raw/genomes/viral/all_viral.fasta
else
    echo "No .fna files to combine"
fi
