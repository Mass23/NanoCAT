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
   conda env create -f NanoCAT.yml  
   conda activate NanoCAT
   ```
2. Process the ONT reads using porechop, chopper, vsearch, and qiime2 using the script: `process_16S.py`
3. Assign taxonomy using the CAT appraoch using the script: `CAT_taxonomies.py`

## Output

After running NanoCAT, the pipeline produces several key files that summarize clustering, classification, and CAT-level aggregation.

### VSEARCH OTU Clustering Output

VSEARCH generates an OTU mapping file: `<folder>/vsearch/otu_clusters.uc`


This file contains:

- ASV-to-OTU cluster assignments  
- The representative centroid (“S” record) sequence for each OTU  
- Dereplicated sequence identifiers and sizes  

This file is used internally by NanoCAT to map read-level taxonomy back onto OTUs.

---

### Aggregated Taxonomy Output

NanoCAT generates one taxonomy table per run, stored at: `<folder>/exports/aggregated_taxonomy_<name>_<mode>.tsv`


Where:

- `<name>` is the user-provided label (`-n`)
- `<mode>` is the taxonomy aggregation method:
  - `centroid` — taxonomy assigned to the OTU representative  
  - `mcCAT` — taxonomy of the read with highest classifier confidence  
  - `hyCAT` — centroid taxonomy replaced only when a more confident read exceeds a multiplier (`-c`)

Each file contains:

| Column        | Description                                                        |
|---------------|--------------------------------------------------------------------|
| Cluster       | OTU identifier (from centroid sequence ID)                         |
| Taxonomy      | Final taxonomy assigned to the OTU                                 |
| Confidence    | Confidence score associated with the assignment                    |

---

### hyCAT Behavior

When using `hyCAT`, NanoCAT applies:

1. Start with the centroid taxonomy  
2. Replace it with the mcCAT taxonomy **only if**:  
   - `Confidence_mcCAT > Confidence_centroid × coefficient`, and  
   - The mcCAT taxonomy is not `"Unassigned"`

This provides a conservative but adaptive alternative to centroid-only classification.


