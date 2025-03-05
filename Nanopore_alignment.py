
import os
import argparse
import pandas as pd
import csv
import subprocess
from pathlib import Path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter and demultiplex reads.")
    parser.add_argument("input", type=str, help="Path to the folder storing the input bam files")
    parser.add_argument("output", type=str, help="Folder path to store the output")
    parser.add_argument("ref_path", type=str, help="Path to the reference sequence fasta")
    #parser.add_argument("ref_path", type=str, help="Path to the reference sequence fasta")

    args = parser.parse_args()

    # Assign to shorter variable names
    input_folder = Path(args.input)
    output_folder = Path(args.output)
    ref_path = Path(args.ref_path)

    if not input_folder.is_dir():
        raise FileNotFoundError(f"Input folder '{input_folder}' does not exist!")
    if not ref_path.is_file():
        raise FileNotFoundError(f"Reference file '{ref_path}' does not exist!")
    output_folder.mkdir(parents=True, exist_ok=True)


    print(f"aligning {input_folder} files to {ref_path}...")
    input_files =  [f.name.replace(".fastq.gz","") for f in input_folder.glob("*.fastq.gz")]
    

    for input_file in input_files:
        output_file = output_folder / input_file  # Use stem of Path object

        print(f"Aligning {input_file}.fastq.gz -> {output_file}.bam")

        #command = f"minimap2 -ax map-ont {ref_path} {input_folder}/{input_file}.fastq.gz > {output_file}.sam" ##Map short accurate genomic reads
        command = f"minimap2 --MD -ax map-ont {ref_path} {input_folder}/{input_file}.fastq.gz | samtools view -bS | samtools sort -o {output_file}.bam"


        subprocess.run(command, shell=True, check=True)
        subprocess.run(f"samtools index {output_file}.bam", shell=True, check=True)

    print("Processing complete!")
