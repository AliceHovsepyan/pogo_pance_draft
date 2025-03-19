import pysam
import os
import argparse
from Bio import SeqIO
import pandas as pd
import csv
from Nanopore_functions import read_cleaning_

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter and demultiplex reads.")
    parser.add_argument("input", type=str, help="Path to the folder storing the input bam files")
    parser.add_argument("output", type=str, help="Folder path to store the output")
    parser.add_argument("ref_path", type=str, help="Path to the reference sequence fasta")

    args = parser.parse_args()

    # Assign to shorter variable names
    input_folder = args.input
    output_folder = args.output
    ref_path = args.ref_path
    cut_n_bases_from_start = 48 ## cut the first 9 bases of the reads

    if not os.path.exists(input_folder):
        raise FileNotFoundError(f"Input folder '{input_folder}' does not exist!")

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)  # Create output folder if missing

    ### read reference sequence
    if not os.path.exists(ref_path):
                    print(f"Reference file {ref_path} does not exist, please check the reference file path")
                    exit()


    ### set the reference
    ref = str(SeqIO.read(ref_path, "fasta").seq)


    all_reads, indels, all_qualitities = read_cleaning_(input_folder, ref, cut_n_bases_from_start)

    ### save cleaned reads
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    with open(f"{output_folder}/cleaned_reads.csv", "w", newline="") as f:
        writer = csv.writer(f)
        for item in all_reads:
            writer.writerow([item]) 

    print(f"Saved cleaned reads to {output_folder}/cleaned_reads.csv")

    with open(f"{output_folder}/ref.csv", "w", newline="") as f:
        writer = csv.writer(f)    
        writer.writerow([ref]) 

    print(f"Saved reference sequence to {output_folder}/ref.csv")

    indels.to_csv(f"{output_folder}/indels.csv") 

    print(f"Saved indels to {output_folder}/indels.csv")
    with open(f"{output_folder}/cleaned_reads_base_qualitities.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(all_qualitities)

    print("Done!")

