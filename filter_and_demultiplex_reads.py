import os
import json
import sys
import subprocess
from pathlib import Path
# from Bio.SeqIO import QualityIO
# import numpy as np
# from matplotlib import pyplot as plt
# import matplotlib.cm as cm
# import gzip
# import glob
# import re
# from DMS_utils import dna_rev_comp, translate_dna2aa
# import pandas as pd
# import seaborn as sns
# import pickle as pkl
# import matplotlib.colors as mcolors
# from scipy import stats
# import os.path
# from matplotlib.lines import Line2D
# import json
# import shutil
# #from evaluation_functions import *
# from functions_ import *
# from plotting import *
# from Bio import SeqIO
# import matplotlib.patches as patches
# from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# from characterization_from_blast_alignments import *

from functions_ import *
from plotting import *


###########define the necessary variables ##########

base_dir = sys.argv[1] if len(sys.argv) > 1 else os.getcwd()

with open(f"{base_dir}/config.json", "r") as file:
    config = json.load(file)

catch_left = config["catch_left"]
catch_right = config["catch_right"]
remove_read_qualities = config["remove_read_qualities"]
Barcodes = config["Barcodes"]
Primer_seq = config["Primer_seq"] 
Primer_out_of_triplets = config["Primer_out_of_triplets"]
variant = config["variant"]
cut_primer_start = config["cut_primer_start"]
cut_BC_seq = config["cut_BC_seq"]
used_Barcodes = config["used_Barcodes"]
Sections = config["Sections"]
amplicon = config["amplicon"]   


quality_score = {
  '!':0, '"':1, '#':2, '$':3, '%':4, '&':5, "'":6, '(':7, ')':8, '*':9,
  '+':10, ',':11, '-':12, '.':13, '/':14, '0':15, '1':16, '2':17, '3':18, '4':19,
  '5':20, '6':21, '7':22, '8':23, '9':24, ':':25, ';':26, '<':27, '=':28, '>':29,
  '?':30, '@':31, 'A':32, 'B':33, 'C':34, 'D':35, 'E':36, 'F':37, 'G':38, 'H':39, 'I':40
}


############ read and demultiplex sequences ###############

## read the sequences, thereby reads are quality filtered, by aborting the read at the first base with a quality score that is in remove_read_qualities list
a_seq, b_seq, _, _, a_ids, b_ids = read_sequences(variant = variant, arbitrary_cutoff_a = False, arbitrary_cutoff_b = False, catch_left=catch_left, catch_right=catch_right, return_qualities_ids=True, quality_score=remove_read_qualities, base_dir = base_dir)

## demultiplex the reads based on the barcodes and the primer sequence
all_reads, all_ids = demultiplex_reads(a_seq, b_seq, ref_gene = None ,Barcodes=Barcodes, Primer_seq=Primer_seq, used_Barcodes = used_Barcodes, Sections = Sections, max_mismatch_primerseq = 5, filter_for_n_mut = False, n_mut_treshold = None, a_ids=a_ids, b_ids=b_ids,  read_len_treshold= None, Primer_out_of_triplets= Primer_out_of_triplets, cut_primer_start=True, cut_BC_seq=True)



############# save as fasta files ###############

## save the demultiplexed reads (R2 reads are saved as reverse complemented reads)
Path(f"{base_dir}/preprocessed/").mkdir(parents = True, exist_ok=True)
Path(f"{base_dir}/references/").mkdir(parents = True, exist_ok=True)


for Bc in used_Barcodes: 

    for section in Sections:

        ref = find_reference_seq(ref_gene=amplicon, Primer_seq=Primer_seq, Section=section, Primer_out_of_triplets=Primer_out_of_triplets) 
        ref_sequences = [SeqRecord(Seq(ref), id = f"{variant}_{section}_ref", description = f"{variant} {section} DNA sequence")]

        for Read_dir in ["R1", "R2"]:

            seqs = all_reads[f"{Bc}_{section}_{Read_dir}"]
            reads = all_reads[f"{Bc}_{section}_{Read_dir}"] if Read_dir == "R1" else [dna_rev_comp(r) for r in all_reads[f"{Bc}_{section}_{Read_dir}"]]

            output_file = f"{base_dir}/preprocessed/{variant}_{Bc}_{section}_Nt_filt_{Read_dir}.fasta"
            sequences = [SeqIO.SeqRecord(Seq(read), id = all_ids[f"{Bc}_{section}_{Read_dir}"][i], description = f"{variant} {Bc} DNA sequence") for i, read in enumerate(reads)]

            count = SeqIO.write(sequences, output_file, "fasta")
            with open(output_file, "w") as output_handle:
                SeqIO.write(sequences, output_handle, "fasta")

            with open(f"{base_dir}/references/{variant}_{Bc}_{section}_Nt_filt_ref.fasta", "w") as output_handle:
                SeqIO.write(ref_sequences, output_handle, "fasta")

            print("Saved %i records to %s" % (count, output_file))
            print(f"Saved reference sequence for {variant} {section} to {base_dir}/preprocessed/{variant}_{Bc}_{section}_Nt_filt_ref.fasta")

    
    


## save ref fastq
# for sec in Sections: 
#     ref = find_reference_seq(ref_gene=amplicon, Primer_seq=Primer_seq, Section=sec, Primer_out_of_triplets=Primer_out_of_triplets) 
#     sequences = [SeqRecord(Seq(ref), id = f"{variant}_{sec}_ref", description = f"{variant} {sec} DNA sequence")]

#     with open(f"{base_dir}/{variant}_{sec}_Nt_ref.fasta", "w") as output_handle:
#         SeqIO.write(sequences, output_handle, "fasta")

#     print(f"Saved reference sequence for {variant} {sec} to {base_dir}/{variant}_{sec}_Nt_ref.fasta")



################ run blast ####################


# Define paths
input_dir =  Path(f"{base_dir}/preprocessed/")      # Directory containing query sequences
reference_dir = Path(f"{base_dir}/references/")   # Directory containing reference sequences
output_dir =  Path(f"{base_dir}/blast/alignments/")   # Directory to store BLAST results
blast_db_dir = Path(f"{base_dir}/blast/db")    # Directory to store BLAST databases

# Ensure output directories exist
output_dir.mkdir(parents=True, exist_ok=True)
blast_db_dir.mkdir(parents=True, exist_ok=True)

# Step 1: Create BLAST databases for each reference sequence
for read_file in input_dir.glob("*.fasta"):
    read_basename = read_file.stem  # Get filename without extension
    db_path = blast_db_dir / read_basename

    print(f"Creating BLAST database for {read_file}...")
    subprocess.run([
        "makeblastdb",
        "-in", str(read_file),
        "-dbtype", "nucl",  # Change to "prot" for protein sequences
        "-out", str(db_path)
    ], check=True)

# Step 2: Run BLAST for each dataset file against their respective reference
for read_file in input_dir.glob("*.fasta"):
    read_basename = read_file.stem

    ref_file = reference_dir / f"{read_basename[:-2]}ref.fasta"

    read_basename = read_file.stem
    db_path = blast_db_dir / read_basename
    output_file = output_dir / f"{read_basename}.out"

    print(f"Running BLAST: aligning reads from {read_file} to {ref_file}...")
    subprocess.run([
        "blastn",  # Use "blastp" for proteins
        "-query", str(ref_file),
        "-db", str(db_path),
        "-out", str(output_file),
        "-outfmt", "15", # Output format 15 is JSON
        "-max_target_seqs", "100000",
    ], check=True)

print("BLAST pipeline completed!")
