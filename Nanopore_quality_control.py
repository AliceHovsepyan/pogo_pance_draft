
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
    #ref_path = args.ref_path
# Step 2: Run BLAST for each dataset file against their respective reference


    print(f"running Nanoplot on {input_folder}...")
    files =  [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith('.bam')]


    subprocess.run([
        "NanoPlot",  # Use "blastp" for proteins
        "-t", "2",
        "--bam", *files,
        "-o", str(output_folder),
        "-f", "pdf",
    ], check=True)
