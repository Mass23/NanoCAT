# 🧬 NanoCAT — Confidence-score based Assignment of Taxonomies

**Version:** v0.1  
**Author:** Massimo Bourquin, November 2025  

---

## What’s NanoCAT?  

Long-read 16S rRNA sequencing with Oxford Nanopore Technologies provides valuable resolution, but elevated error rates can complicate taxonomic classification and reduce assignment accuracy. NanoCAT addresses this challenge by integrating read-level classifier confidence, enabling more reliable OTU-level taxonomy that better reflects the underlying data.

This is the ONT 16S pipeline of the MACE laboratory (EPFL). 

---

## What’s Inside  

- Scripts for the full ONT 16S pipeline: from read preprocessing to taxonomy assignment.
- The CAT approach, including:
  
--> **mcCAT** — Selects the taxonomy of the **single read with the highest confidence** in each OTU.  
--> **hyCAT** — Uses the highest-confidence read *only if* it beats what centroid-based classification would have called.
  
- Support files to reproduce the analyses from our study ("Beyond the centroid..."): conda environment specs (`beyond_the_centroid.yaml`) and R analysis scripts (`beyond_the_centroid_analyses.R`).

---

## Installation & Usage

1. Create the conda environment:  
   ```bash
   conda env create -f CAT.yml  
   conda activate NanoCAT
   ```
2. Process the ONT reads using porechop, chopper, vsearch, and qiime2 using the script: `process_16S.py`
3. Assign taxonomy using the CAT appraoch using the script: `CAT_taxonomies.py`

## Output
