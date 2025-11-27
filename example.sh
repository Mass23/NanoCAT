python3 process_16S_nanopore.py \
    -f test/raw_reads \
    -n test/test \
    -m test/metadata.tsv \
    -t 2 \
    -c test/classifier/silva-138-99-nb-classifier.qza

# CAT using hyCAT
python3 CAT_taxonomy.py -f test/test_results/ -t test/test_results/exports/taxonomy.tsv -n test_silva_leale2025 -m hyCAT -c 1.2

# CAT using mcCAT
python3 CAT_taxonomy.py -f test/test_results/ -t test/test_results/exports/taxonomy.tsv -n test_silva_leale2025 -m mcCAT
