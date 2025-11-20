# NanoCAT 
#### v0.1, Author: Massimo Bourquin, October 2025

## Description
This is the ONT 16S pipeline of the MACE laboratory (EPFL). It contains the code for the CAT (Confidence score based Assignment of Taxonomies), and other scripts to generate the necessary data.

## Installation & Usage

1. Install the conda environment in: `CAT.yml`
2. Process your data using porechop, chopper, vsearch, and qiime2 using the script: `process_16S.py`
3. Assign taxonomy using the CAT appraoch using the script: `CAT_taxonomies.py`

to do: upload scripts and describe output

## Others

The `beyond_the_centroid.yaml` and `beyond_the_centroid_analyses.R` are the conda environment and the R script that can be used to replicate the study: "Beyond the centroid: enhancing taxonomic profiling of Nanopore 16S rRNA data using classification confidence scores".

## Output
