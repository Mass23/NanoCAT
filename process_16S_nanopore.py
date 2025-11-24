import os
import argparse
import subprocess
import pandas as pd # type: ignore
import datetime
import glob

################################################################################
#################           FUNCTIONS           ################################
################################################################################
def add_to_log(results_folder_name, text):
    with open(f'{results_folder_name}/log.txt', 'w') as log:
        log.write(text + '\n\n')

# Preprocessing part
def create_result_folder(results_folder_name):
    """
    Creates a folder results, if the folder already exists, does not overwrites it.
    """
    if not os.path.exists(results_folder_name):
        os.makedirs(results_folder_name)

    add_to_log(results_folder_name, f"Log file for the run {results_folder_name}, time and date: {datetime.datetime.now().strftime('%I:%M%p on %B %d, %Y')}" + '\n\n')

def print_env_summary(results_folder_name):
    """
    Print the conda environment list (packages and softwares with version numbers).
    """
    subprocess.call(f'conda list > {results_folder_name}/list_conda.txt', shell = True)
    add_to_log(results_folder_name, 'Creating the list_conda.txt file that summarises the conda env.' + '\n\n')

def load_metadata(metadata_path, extension):
    """
    Load the metadata as a pandas data frame.
    """
    print(extension)
    if extension == 'tsv':
        metadata = pd.read_csv(metadata_path, sep='\t', header=0)
    elif extension == 'csv':
        metadata = pd.read_csv(metadata_path, sep=',', header=0)
        print(metadata.head())
    else:
        print('Only tsv or csv files are accepted!!!')
    return(metadata)

def list_fastq_files(folder_path):
    """
    Lists all folders in the given folder (= samples fastq files).
    """
    try:
        # Check if the path is a directory
        if not os.path.isdir(folder_path):
            print(f"The path '{folder_path}' is not a directory.")
            return

        # List all entries in the directory
        fastq_files = glob.glob(f'{folder_path}/*.fastq.gz')

        # Filter out files, keep only directories
        filtered_fastq_files = [fastq for fastq in fastq_files if fastq != 'unclassified']
        return(filtered_fastq_files)

    except Exception as e:
        print(f"An error occurred: {e}")

def list_files(folder_path):
    """
    List the fastq.gz files in the specifiec raw reads directory.
    """
    files = glob.glob(f'{folder_path}/*.fastq.gz')
    samples = [file.replace('.fastq.gz','') for file in files]
    samples = [os.path.basename(file) for file in samples]
    samples = [file.split('_')[-1] for file in samples]
    filtered_samples = [sample for sample in samples if not 'unclassified' in sample]
    return(filtered_samples)

def list_samples(output_path):
    files = glob.glob(f'{output_path}/raw_data/*_chopped.fastq.gz')
    samples = [os.path.basename(file).replace('_chopped.fastq.gz','') for file in files]
    return(samples)

def check_metadata_samples(metadata, samples, results_folder_name):
    """
    Check that the metadata agrees with the sample names listed.
    """
    metadata_samples = set(metadata['Barcode']) # samples in the metadata
    files_samples = set(samples) # samples in the reads' files
    print(files_samples)

    print('There are %d samples in the metadata file'%(len(metadata_samples)))
    print("There are %d samples in the reads' files"%(len(files_samples)))
    print("Barcode in the metadata but not in the reads' files:")
    print(metadata_samples.difference(files_samples))

    with open(f'{results_folder_name}/log.txt', 'a') as log:
        log.write('There are %d samples in the metadata file'%(len(metadata_samples)) + '\n\n')
        log.write("There are %d samples in the reads' files"%(len(files_samples)) + '\n\n')
        log.write("Barcode in the metadata but not in the reads' files:" + '\n\n')
        log.write(str(metadata_samples.difference(files_samples)) + '\n\n')

    out_list = [sample for sample in list(metadata_samples) if sample in files_samples]
    return metadata.loc[metadata['Barcode'].isin(out_list)], out_list

def concatenate_files(folder_path, metadata, samples, results_folder_name):
    """
    Takes folder paths, metadata files and samples list and concatenate the data
    for each sample in the list, + renames it to the sample name using metadata.
    """
    os.makedirs('%s/raw_data' % (results_folder_name))
    new_samples = []
    for sample in samples:
        new_sample = metadata.loc[metadata['Barcode'] == sample, 'Sample'].values[0]
        new_path = f'{results_folder_name}/raw_data/{new_sample}.fastq.gz'
        args = f'cat {folder_path}/*{sample}.fastq.gz > {new_path}'
        subprocess.call(args, shell = True)
        new_samples.append(new_sample)
    return(new_samples)

# Reads preprocessing part
def run_porechop(samples, threads, results_folder_name):
    add_to_log(results_folder_name, f'Running porechop...')

    for sample in samples:
        reads_in = f'{results_folder_name}/raw_data/{sample}.fastq.gz'
        reads_out = f'{results_folder_name}/raw_data/{sample}_porechopped.fastq.gz'
        args = f'porechop --threads {str(threads)} -i {reads_in} -o {reads_out}'
        subprocess.call(args, shell = True)
        add_to_log(results_folder_name, args)

def run_chopper(samples, threads, results_folder_name, qual_threshold):
    add_to_log(results_folder_name, f'Running chopper...')

    for sample in samples:
        reads_in = f'{results_folder_name}/raw_data/{sample}_porechopped.fastq.gz'
        reads_out = f'{results_folder_name}/raw_data/{sample}_chopped.fastq.gz'
        args = f'gunzip -c {reads_in} | chopper -q {str(qual_threshold)} --maxlength 1800 --minlength 1200 --threads {str(threads)} | gzip > {reads_out}'
        subprocess.call(args, shell = True)
        add_to_log(results_folder_name, args)

# vsearch part
def run_vsearch(results_folder_name, samples, threads, perc_identity):
    os.makedirs(f'{results_folder_name}/vsearch')
    os.makedirs(f'{results_folder_name}/vsearch/drep_data')

    for sample in samples:
        fastq_in = f'{results_folder_name}/raw_data/{sample}_chopped.fastq.gz'
        fasta_out = f'{results_folder_name}/vsearch/drep_data/{sample}_raw.fasta'
        copy_args = f'seqtk seq -a {fastq_in} > {fasta_out}'
        subprocess.call(copy_args, shell = True)
        add_to_log(results_folder_name, copy_args)

        args_1 = f'vsearch --derep_fulllength {fasta_out} --strand both --output {results_folder_name}/vsearch/drep_data/{sample}_derep.fasta --sizeout --uc {results_folder_name}/vsearch/drep_data/{sample}s.derep.uc --relabel {sample}.s --fasta_width 0'
        subprocess.call(args_1, shell=True)
        add_to_log(results_folder_name, args_1)

    args_2 = f'cat {results_folder_name}/vsearch/drep_data/*derep.fasta > {results_folder_name}/vsearch/drep_data/all_derep.fasta'
    subprocess.call(args_2, shell=True)
    add_to_log(results_folder_name, args_2)

    #args_3 = f'vsearch --derep_fulllength {results_folder_name}/vsearch/drep_data/all_derep.fasta --sizein --sizeout --fasta_width 0 --uc {results_folder_name}/vsearch/merged.derep.uc --output {results_folder_name}/vsearch/merged.derep.fasta'
    #subprocess.call(args_3, shell = True)
    #add_to_log(results_folder_name, args_3)

    args_4 = f'vsearch --cluster_size {results_folder_name}/vsearch/drep_data/all_derep.fasta --threads {str(threads)} --id {str(perc_identity)} --strand both --sizein --sizeout --fasta_width 0 --uc {results_folder_name}/vsearch/otu_clusters.uc --centroids {results_folder_name}/vsearch/otu_centroids.fasta  --otutabout {results_folder_name}/vsearch/otu_table.tsv'
    subprocess.call(args_4, shell = True)
    add_to_log(results_folder_name, args_4)

# Qiime 2 part
#def create_manifest(results_folder_name):
#    os.makedirs(f'{results_folder_name}/qiime2')
#    full_path = os.getcwd()
#    manifest = pd.DataFrame(columns=['sample-id', 'absolute-filepath'])
#    
#    sample_files = glob.glob(f'{results_folder_name}/raw_data/*_chopped.fastq.gz')
#    samples = [sample.split('/')[2].replace('_chopped.fastq.gz','') for sample in sample_files]
#    print(samples)
#    for sample in samples:
#        manifest.loc[len(manifest.index)] = [sample, f'{full_path}/{results_folder_name}/raw_data/{sample}_chopped.fastq.gz'] 
#    manifest.to_csv(f'{results_folder_name}/qiime2/qiime2_manifest.tsv', sep='\t', index=False) 
#
#    with open(f'{results_folder_name}/log.txt', 'a') as log:
#        log.write('Qiime2 manifest created...' + '\n\n')

# taxonomy part
def taxonomy_qiime2(results_folder_name, classifier_path, threads):
    os.makedirs(f'{results_folder_name}/qiime2')

    args_1 = f"qiime tools import --input-path {results_folder_name}/vsearch/drep_data/all_derep.fasta --output-path {results_folder_name}/qiime2/sequences.qza --type 'FeatureData[Sequence]'"
    subprocess.call(args_1, shell = True)
    add_to_log(results_folder_name, args_1)

    if ',' in classifier_path:
        classifiers = classifier_path.split(',')
        for classifier in classifiers:
            print(classifier)
            classifier_name = classifier.split('/')[-1].replace('.qza','')
            args_2 = f'qiime feature-classifier classify-sklearn --p-n-jobs {threads} --i-reads {results_folder_name}/qiime2/sequences.qza --i-classifier {classifier} --o-classification {results_folder_name}/qiime2/taxonomy-classification-{classifier_name}.qza'
            subprocess.call(args_2, shell = True)
            add_to_log(results_folder_name, args_2)

            args_3 = f'qiime tools export --input-path {results_folder_name}/qiime2/taxonomy-classification-{classifier_name}.qza --output-path {results_folder_name}/exports/taxonomy-classification-{classifier_name}'
            subprocess.call(args_3, shell = True)
            add_to_log(results_folder_name, args_3)
    else:
        args_2 = f'qiime feature-classifier classify-sklearn --p-n-jobs {threads} --i-reads {results_folder_name}/qiime2/sequences.qza --i-classifier {classifier_path} --o-classification {results_folder_name}/qiime2/taxonomy-classification.qza'
        subprocess.call(args_2, shell = True)
        add_to_log(results_folder_name, args_2)

        args_3 = f'qiime tools export --input-path {results_folder_name}/qiime2/taxonomy-classification.qza --output-path {results_folder_name}/exports'
        subprocess.call(args_3, shell = True)
        add_to_log(results_folder_name, args_3)

################################################################################
#################             MAIN             #################################
################################################################################
def main():
    # Create an argument parser
    parser = argparse.ArgumentParser()

    # Add the folder path argument
    parser.add_argument("-f", "--folder", type=str,
                        help="Path to the reads' folder.", required=True)
    parser.add_argument("-n", "--name", type=str,
                        help="Name of the results folder (_results will be added at the end), please keep the output in the processed_data folder.", required=True)
    parser.add_argument("-m", "--metadata_file", type=str,
                        help="Path to the metadata file, if it is csv or tsv, the parser recognises the separator based on the extension.", required=True)
    parser.add_argument("-t", "--threads", type=str,
                        help="Number of threads to use for multiprocessing-compatible tasks.", required=True)
    parser.add_argument("-c", "--classifier", type=str,
                        help="Path of the classifier to use, if you want to use several, provide a comma separated list of paths.", required=True)
    parser.add_argument("--skippreprocessing", action='store_true',
                        help="To add if you want to skip preprocessing.")
    parser.add_argument("--onlypreprocessing", action='store_true',
                        help="To add if you want to skip the vsearch and qiime2 parts (only preprocessing).")
    parser.add_argument("--perc_identity", type=int, required=False, default='97',
                        help="To add if you want to perform OTU clustering at a different threshold, default is 97%.")
    parser.add_argument("--qual_threshold", type=int, required=False, default='12',
                        help="To add if you want to change the quality threshold for chopper, default is 12.")


    # Parse arguments
    args = parser.parse_args()
    out_folder = f'{args.name}_results'

    if args.skippreprocessing is False:
        # Create results folder, print the environment summary, load the metadata
        # and list the samples to process
        create_result_folder(out_folder)
        print_env_summary(out_folder)
        ext = args.metadata_file.split('.')[-1]
        metadata = load_metadata(args.metadata_file, ext)
        samples = list_files(args.folder)
        metadata, samples_to_process = check_metadata_samples(metadata, samples, out_folder)

        # Concatenate files belonging to the same sample in the new directory,
        # run porechop and chopper
        samples_names = concatenate_files(args.folder, metadata, samples_to_process, out_folder)
        run_porechop(samples_names, args.threads, out_folder)
        run_chopper(samples_names, args.threads, out_folder, args.qual_threshold)
    else:
        samples_names = list_samples(out_folder)
        print('Samples in raw_data folder:')
        print(samples_names)


    if args.onlypreprocessing is False:
        # Create the Qiime manifest, run qiime analysis
        run_vsearch(out_folder, samples_names, args.threads, args.perc_identity/100)
        #import_qiime2(out_folder)
        taxonomy_qiime2(out_folder, args.classifier, args.threads)

if __name__ == "__main__":
    main()
