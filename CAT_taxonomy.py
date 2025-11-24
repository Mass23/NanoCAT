import numpy as np
import pandas as pd
import argparse
import os
import random
random.seed(10)

###############################################################################
# ---------------------------  TAXONOMY HELPERS  ------------------------------
###############################################################################

def GetCentroidTaxonomy(tab):
    """Return taxonomy assigned to the OTU representative (S record)."""
    return tab.loc[tab.Type == 'S'][['Cluster', 'Taxonomy', 'Confidence']]


def GetMcCATTaxonomy(tab):
    """Return mcCAT taxonomy = taxonomy of the read with highest confidence."""
    idx = tab.groupby("Cluster").Confidence.idxmax()
    return tab.loc[idx][['Cluster', 'Taxonomy', 'Confidence']]


def GetHyCATTaxonomy(tab, coeff):
    """
    hyCAT:
      - Start from centroid taxonomy
      - If the most-confident read exceeds (centroid_conf * coeff),
        replace taxonomy.
    """
    # Extract centroid-level taxonomy
    centroid = GetCentroidTaxonomy(tab).rename(columns={
        "Taxonomy": "Taxonomy.centroid",
        "Confidence": "Confidence.centroid"
    })

    # Extract mcCAT taxonomy
    mc = GetMcCATTaxonomy(tab).rename(columns={
        "Taxonomy": "Taxonomy.mcCAT",
        "Confidence": "Confidence.mcCAT"
    })

    # Merge
    merged = pd.merge(centroid, mc, on="Cluster", how="left")

    # Default hyCAT = centroid
    merged["Taxonomy.hyCAT"] = merged["Taxonomy.centroid"]
    merged["Confidence.hyCAT"] = merged["Confidence.centroid"]

    # Apply replacement rule
    mask = (
        (merged["Confidence.mcCAT"] > (merged["Confidence.centroid"] * coeff)) &
        (merged["Taxonomy.mcCAT"] != "assigned")
    )

    merged.loc[mask, "Taxonomy.hyCAT"] = merged.loc[mask, "Taxonomy.mcCAT"]
    merged.loc[mask, "Confidence.hyCAT"] = merged.loc[mask, "Confidence.mcCAT"]

    return merged[['Cluster', 'Taxonomy.hyCAT', 'Confidence.hyCAT']].rename(
        columns={"Taxonomy.hyCAT": "Taxonomy", "Confidence.hyCAT": "Confidence"}
    )


###############################################################################
# ------------------------------  PARSING  ------------------------------------
###############################################################################

def ParseTables(otufile, taxfile):
    otutab = pd.read_csv(
        otufile, sep='\t', header=None,
        usecols=[0, 1, 8, 9],
        names=['Type', 'Cluster', 'ASVsize', 'Centroid']
    )
    taxtab = pd.read_csv(
        taxfile, sep='\t', header=0,
        names=['ASVsize', 'Taxonomy', 'Confidence']
    )

    # Extract centroid ID
    otutab['Centroid'] = otutab.index.map(
        lambda x: otutab.loc[x, 'Centroid'].split(';')[0]
        if otutab.loc[x, 'Centroid'] != "*"
        else otutab.loc[x, 'ASVsize'].split(';')[0]
    )

    # Merge with taxonomy
    fulltab = pd.merge(otutab, taxtab, on='ASVsize')
    fulltab = fulltab[fulltab.Type != 'C']

    # Clean cluster identifiers
    fulltab['Cluster'] = fulltab['Centroid'].apply(lambda x: x.split(';')[0])

    # Extract ASV and abundance
    fulltab[['ASV', 'Count']] = fulltab['ASVsize'].str.split(';size=', n=1, expand=True)
    fulltab['Count'] = fulltab['Count'].astype(int)

    # Compute per-OTU abundance
    fulltab['ClusterCount'] = fulltab['Count'].groupby(fulltab['Cluster']).transform('sum')

    return fulltab


###############################################################################
# ------------------------------- MAIN ----------------------------------------
###############################################################################

def main():
    parser = argparse.ArgumentParser(
        description="Aggregates OTU taxonomy using centroid, mcCAT, or hyCAT methods."
    )

    parser.add_argument("-f", "--folder", required=True, type=str,
                        help="Path to the results folder.")
    parser.add_argument("-t", "--taxonomy", required=True, type=str,
                        help="Taxonomy file returned by classifier.")
    parser.add_argument("-n", "--name", required=True, type=str,
                        help="Name to append to output file.")
    parser.add_argument("-m", "--mode", default='centroid', type=str,
                        help="Aggregation mode: centroid / mcCAT / hyCAT.")
    parser.add_argument("-c", "--coeff", default=1.0, type=float,
                        help="Multiplier threshold for hyCAT (default=1.0).")

    args = parser.parse_args()

    fulltab = ParseTables(f"{args.folder}vsearch/otu_clusters.uc", args.taxonomy)

    # --------- SELECT AGGREGATION METHOD ---------
    if args.mode == 'centroid':
        agg_tax = GetCentroidTaxonomy(fulltab)

    elif args.mode == 'mcCAT':
        agg_tax = GetMcCATTaxonomy(fulltab)

    elif args.mode == 'hyCAT':
        agg_tax = GetHyCATTaxonomy(fulltab, coeff=args.coeff)

    else:
        raise ValueError("Mode not recognized. Use: centroid / mcCAT / hyCAT")

    # --------- OUTPUT ---------
    outfile = f"{args.folder}exports/aggregated_taxonomy_{args.name}_{args.mode}.tsv"
    agg_tax.to_csv(outfile, sep="\t", index=False, encoding="utf-8")
    print(f"Saved: {outfile}")


if __name__ == "__main__":
    main()
