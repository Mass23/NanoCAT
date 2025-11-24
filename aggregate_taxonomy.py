import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import random
import argparse
import subprocess
random.seed(10)

def GetCentroidTaxonomy(tab):
    return(tab.loc[tab.Type == 'S'][['Cluster', 'Taxonomy', 'Confidence']])

def GetMostConfidentTaxonomy(tab):
    return(tab.loc[tab.groupby(["Cluster"]).Confidence.idxmax()][['Cluster', 'Taxonomy', 'Confidence']])

def AssignedConsensusTaxonomy(taxonomies):
    assigned = taxonomies[taxonomies != 'Unassigned']
    if len(assigned) == 0:
        return('Unassigned')
    else:
        alltax = assigned.str.split('; ').explode(ignore_index=True)
        doma = alltax[alltax.str.startswith('d__')].mode().to_string(index=False)
        phyl = alltax[alltax.str.startswith('p__')].mode().to_string(index=False)
        clas = alltax[alltax.str.startswith('c__')].mode().to_string(index=False)
        orde = alltax[alltax.str.startswith('o__')].mode().to_string(index=False)
        fami = alltax[alltax.str.startswith('f__')].mode().to_string(index=False)
        genu = alltax[alltax.str.startswith('g__')].mode().to_string(index=False)
        spec = alltax[alltax.str.startswith('s__')].mode().to_string(index=False)
        return('; '.join([doma,phyl,clas,orde,fami,genu,spec]).replace(' Series([], )',''))

def GetAssignedConsensusTaxonomy(tab):
    taxonomies = tab.groupby(["Cluster"]).Taxonomy.apply(AssignedConsensusTaxonomy)
    clusters = tab.groupby(["Cluster"]).Cluster.first()
    return(pd.DataFrame({'Cluster': clusters.tolist(), 'Taxonomy': taxonomies.tolist()}).reset_index())

def AlignCluster(sequences, cpus):
    with open('tmp.fasta', 'w') as handle:
        SeqIO.write(sequences.values(), handle, 'fasta')
    subprocess.call(f'mafft --thread {str(cpus)} tmp.fasta > tmp_aln.fasta', shell=True)
    align = AlignIO.read('tmp_aln.fasta', "fasta")
    os.remove('tmp.fasta')
    os.remove('tmp_aln.fasta')
    return(align)

def GetConsensusAlignmentTaxonomy(tab, fasta_file, cpus):
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    sequences = {key.split(';')[0]: value for key, value in sequences.items()}
    clusters = tab.Cluster.unique().tolist()
    
    consensus_fasta = []
    for cluster in clusters:
        print(cluster)
        to_keep = tab[tab.Cluster == cluster]['ASV'].tolist()
        if len(to_keep) > 50:
            to_keep = random.sample(to_keep, 50)
        cluster_sequences = {k: sequences[k] for k in to_keep}
        print(len(cluster_sequences))
        alignment = AlignCluster(cluster_sequences, cpus)
        summary_align = AlignInfo.SummaryInfo(alignment)
        consensus = summary_align.dumb_consensus(threshold=0.5)
        consensus_fasta.append(SeqRecord(Seq(consensus), str(cluster)))
    SeqIO.write(consensus_fasta, fasta_file.replace('all_derep','all_derep_consensus'), "fasta")



def ParseTables(otufile, taxfile):
    otutab = pd.read_csv(otufile, sep = '\t', header = None, usecols = [0, 1, 8, 9], names = ['Type', 'Cluster', 'ASVsize', 'Centroid'])
    taxtab = pd.read_csv(taxfile, sep='\t', header=0, names=['ASVsize', 'Taxonomy', 'Confidence'])

    otutab['Centroid'] = otutab.index.map(lambda x: otutab.loc[x, 'Centroid'].split(';')[0] if otutab.loc[x, 'Centroid'] != "*" else otutab.loc[x, 'ASVsize'].split(';')[0])
    fulltab = pd.merge(otutab, taxtab, on='ASVsize')
    fulltab = fulltab[fulltab.Type != 'C']

    fulltab['Cluster'] = fulltab['Centroid'].apply(lambda x: x.split(';')[0]) #'clu_' + fulltab['Cluster'].astype(str)
    fulltab[['ASV', 'Count']] = fulltab['ASVsize'].str.split(';size=', n=1, expand=True)
    fulltab.Count = fulltab.Count.astype('int64')
    fulltab['ClusterCount'] = fulltab['Count'].groupby(fulltab['Cluster']).transform('sum')
    #fulltab = fulltab[fulltab.ClusterCount > 1]

    return(fulltab)

def main():
    # Create an argument parser
    parser = argparse.ArgumentParser(description="Aggregates OTU taxonomy from ASV taxonomy classification.")

    # Add the folder path argument
    parser.add_argument("-f", "--folder", type=str, help="Path to the results folder as a string.", required=True)
    parser.add_argument("-t", "--taxonomy", type=str, help="Name of the taxonomy file to process.", required=True)
    parser.add_argument("-n", "--name", type=str, help="Name to append to the output file.", required=True)
    parser.add_argument("-m", "--mode", type=str, default='centroid', help="What method to use to aggregate taxonomy: centroid/mostconfident/taxonconsensus/alignconsensus, default is centroid.")              
    parser.add_argument("-c", "--cpus", type=int, default=1, help="Number of CPUs to use.")      
    

    # Parse arguments
    args = parser.parse_args()

    fulltab = ParseTables(f'{args.folder}vsearch/otu_clusters.uc', args.taxonomy)

    if args.mode == 'centroid':
        agg_tax = GetCentroidTaxonomy(fulltab)
    
    elif args.mode == 'mostconfident':
        agg_tax = GetMostConfidentTaxonomy(fulltab)
    
    elif args.mode == 'taxonconsensus':
        agg_tax = GetAssignedConsensusTaxonomy(fulltab)

    elif args.mode == 'alignconsensus':
        agg_tax = GetConsensusAlignmentTaxonomy(fulltab, f'{args.folder}vsearch/drep_data/all_derep.fasta', args.cpus)
    
    else:
        print('Mode not recognised, please choose one of centroid, mostconfident, taxonconsensus & alignconsensus')

    if args.mode == 'alignconsensus':
        return

    else:
        agg_tax.to_csv(f'{args.folder}exports/aggregated_taxonomy_{args.name}_{args.mode}.tsv', encoding='utf-8', index=False)

if __name__ == "__main__":
    main()
