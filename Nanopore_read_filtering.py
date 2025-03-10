
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
    input_files =  [f for f in os.listdir(input_folder) if f.endswith(".fastq.gz")]

    
    # Define input and output directories

    # Loop through all `.fastq.gz` files in the input directory
    for input_file in  input_files:
        print(f"Processing {input_file}")
        
        output_file = f"{output_folder}/{input_file}"

        # Construct and run the command
        command = f"chopper -q 20 --minlength 1800 --maxlength 2100 -i {input_folder}/{input_file} | gzip > {output_file}"
        print(f"Processing {input_file} -> {output_file}")

        subprocess.run(command, shell=True, check=True)

    print("Processing complete!")
