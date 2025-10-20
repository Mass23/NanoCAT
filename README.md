# NanoCAT 
#### v0.1, Author: Massimo Bourquin, October 2025

## Description
This is the ONT 16S pipeline of the MACE laboratory (EPFL). It contains the code for the CAT (Confidence score based Assignment of Taxonomies), and other scripts to generate the necessary data.

## Installation & Usage

1. Install the conda environment in: `env/CAT.yml`
2. Process your data using porechop, chopper, vsearch, and qiime2 using the script: `process_16S.py`
3. Assign taxonomy using the CAT appraoch using the script: `CAT_taxonomies.py`

## Output
