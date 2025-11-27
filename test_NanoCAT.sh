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
SRA_ID="DRR225044"

echo "[INFO] Downloading FASTQ for ${SRA_ID}…"
# The --split-3 ensures single-end/single reads are handled
fastq-dump --split-3 --outdir test/raw_reads "${SRA_ID}"
gzip test/raw_reads/DRR225044.fastq

### -----------------------------
### 3. Create metadata.tsv
### -----------------------------
echo "[INFO] Creating metadata file…"

printf "Sample\tBarcode\n%s\t%s\n" "10 Strains even mix" "$SRA_ID" > test/metadata.tsv

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
    -t 2 \
    -c test/classifier/silva-138-99-nb-classifier.qza

echo "[INFO] Test completed successfully!"
