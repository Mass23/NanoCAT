#!/usr/bin/env bash
set -euo pipefail

### -----------------------------
### 1. Create test folders
### -----------------------------
echo "[INFO] Creating test directories…"
mkdir -p test/raw_reads
mkdir -p test/classifier

### -----------------------------
### 2. Download reads using fastq-dump
### -----------------------------
SRA_ID="SRR28553412"

echo "[INFO] Downloading FASTQ for ${SRA_ID}…"
# The --split-3 ensures single-end/single reads are handled
fastq-dump --split-3 --outdir test/raw_reads "${SRA_ID}"
gzip test/raw_reads/SRR28553412.fastq

### -----------------------------
### 3. Create metadata.tsv
### -----------------------------
echo "[INFO] Creating metadata file…"

printf "Sample\tBarcode\n%s\t%s\n" "sample_from_leale2024" "$SRA_ID" > test/metadata.tsv

echo "[INFO] Metadata created at test/metadata.tsv"

### -----------------------------
### 4. Download classifier
### -----------------------------
CLASSIFIER_URL="https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza"

echo "[INFO] Downloading QIIME2 classifier…"
wget -O test/classifier/silva-138-99-nb-classifier.qza "${CLASSIFIER_URL}"

### -----------------------------
### 5. Activate conda env NanoCAT
### -----------------------------
echo "[INFO] Activating NanoCAT conda environment…"

# Ensure conda is available in the shell
# (this line is needed in non-interactive scripts)
source "$(conda info --base)/etc/profile.d/conda.sh"

conda activate NanoCAT

### -----------------------------
### 6. Run the pipeline
### -----------------------------
echo "[INFO] Running process_16S_nanopore.py…"

python3 process_16S_nanopore.py \
    -f test/raw_reads \
    -n test/test \
    -m test/metadata.tsv \
    -t 2 --perc_identity 90 --qual_threshold 10 \
    -c test/classifier/silva-138-99-nb-classifier.qza

echo "[INFO] Running mcCAT and hyCAT…"
python3 CAT_taxonomy.py -f test/test_results/ -t test/test_results/exports/taxonomy.tsv -n test_silva_leale2025 -m hyCAT -c 1.2
python3 CAT_taxonomy.py -f test/test_results/ -t test/test_results/exports/taxonomy.tsv -n test_silva_leale2025 -m mcCAT

echo "[INFO] Checking MD5 integrity of final output files…"

MD5FILE="test/expected_md5s.txt"
rm -f "$MD5FILE"

# Write expected md5 checksums
printf "%s  %s\n" "3ecccf677205da601bb2e1d7c1dbaff7" "test/test_results/vsearch/otu_table.tsv" >> "$MD5FILE"
printf "%s  %s\n" "a05d8086f484fbd6b387facb6029d811" "test/test_results/exports/aggregated_taxonomy_test_silva_leale2025_hyCAT.tsv" >> "$MD5FILE"
printf "%s  %s\n" "9af98bd9901486d0b1b7548759d051ac" "test/test_results/exports/aggregated_taxonomy_test_silva_leale2025_mcCAT.tsv" >> "$MD5FILE"

echo "[INFO] Expected MD5 file created:"
cat "$MD5FILE"

# Verify MD5 checks — Linux or macOS
if command -v md5sum >/dev/null 2>&1; then
    echo "[INFO] Running md5sum check…"
    md5sum -c "$MD5FILE" || { echo "[ERROR] MD5 verification failed!"; exit 1; }
else
    echo "[INFO] md5sum not found, using macOS md5 fallback."
    while read -r md5 path; do
        computed=$(md5 -r "$path" | awk '{print $1}')
        if [[ "$computed" != "$md5" ]]; then
            echo "[ERROR] MD5 mismatch for $path"
            echo "  Expected: $md5"
            echo "  Got:      $computed"
            exit 1
        else
            echo "[OK] $path"
        fi
    done < "$MD5FILE"
fi

echo "[INFO] All MD5 checks passed."


echo "[INFO] Test completed successfully!"
