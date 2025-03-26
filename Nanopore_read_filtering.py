
import os
import argparse
import pandas as pd
import csv
import subprocess

### use chopper for quality filtering and filtering after read length
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter and demultiplex reads.")
    parser.add_argument("input", type=str, help="Path to the folder storing the input bam files")
    parser.add_argument("output", type=str, help="Folder path to store the output")

    args = parser.parse_args()

    # Assign to shorter variable names
    input_folder = args.input
    output_folder = args.output
    if not os.path.exists(input_folder):
        raise FileNotFoundError(f"Input folder '{input_folder}' does not exist!")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)


    print(f"running Nanoplot on {input_folder}...")
    input_files =  [f for f in os.listdir(input_folder) if f.endswith(".fastq.gz")]

    for input_file in  input_files:
        print(f"Processing {input_file}")
        
        output_file = f"{output_folder}/{input_file}"

        # Construct and run the command
        command = f"chopper -q 20 --minlength 1800 --maxlength 2200  -i {input_folder}/{input_file} | gzip > {output_file}" #--headcrop 50 --endcrop 50 --minlength 1800 --maxlength 2100 
        print(f"Processing {input_file} -> {output_file}")

        subprocess.run(command, shell=True, check=True)

    print("Processing complete!")
