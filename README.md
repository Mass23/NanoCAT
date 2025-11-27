# 🧬 NanoCAT — Confidence-score based Assignment of Taxonomies

**Version:** v0.1  
**Author:** Massimo Bourquin, November 2025  

/For now tested on linux/

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

First, install the environment with the yml file available in this repository.
   ```bash
   conda create --file NanoCAT.yml
   conda activate NanoCAT
   ```

The installation can be tested using the bash script: `test_NanoCAT.sh`. The test will download one SRR sequencing data from the internet, one qiime2 classifier and process that sample using NanoCAT (`process_16S_nanopore.py` and then `CAT_taxonomy.py`), and then test the output files using md5sum checks.

3. Process the ONT reads using porechop, chopper, vsearch, and qiime2 using the script: `process_16S_nanopore.py`
4. Assign taxonomy using the CAT appraoch using the script: `CAT_taxonomy.py`

## Output

After running NanoCAT, the pipeline produces several key files that summarize clustering, classification, and CAT-level aggregation. To analyse the dataset the important ones are:

1. The OTU table in output_name`_results/vsearch/otu_table.tsv`

2. After running the CAT aggregation script, NanoCAT generates one taxonomy table per run, stored at: `<folder>/exports/aggregated_taxonomy_<name>_<mode>.tsv`

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

This provides a conservative but adaptive alternative to centroid-only and most-confident classifications.


