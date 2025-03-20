
import os
import argparse
import pandas as pd
import csv
import subprocess

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter and demultiplex reads.")
    parser.add_argument("input", type=str, help="Path to the folder storing the input bam files")
    parser.add_argument("output", type=str, help="Folder path to store the output")
    #parser.add_argument("ref_path", type=str, help="Path to the reference sequence fasta")

    args = parser.parse_args()

    # Assign to shorter variable names
    input_folder = args.input
    output_folder = args.output
    if not os.path.exists(input_folder):
        raise FileNotFoundError(f"Input folder '{input_folder}' does not exist!")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)


    print(f"running Nanoplot on {input_folder}...")
    filetype = '.fastq' if '.fastq.gz' in os.listdir(input_folder)[0] else '.bam'
    fileend = '.fastq.gz' if filetype == '.fastq' else '.bam'
    files =  [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith(fileend)]#bam fastq.gz


    subprocess.run([
        "NanoPlot",  # Use "blastp" for proteins
        "-t", "2",
        f"--{filetype[1:]}", *files, #bam
        "-o", str(output_folder),
        "-f", "pdf",
    ], check=True)
