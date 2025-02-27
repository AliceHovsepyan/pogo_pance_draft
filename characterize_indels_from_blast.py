import os
from Bio.SeqIO import QualityIO
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import gzip
import glob
import re
from DMS_utils import dna_rev_comp, translate_dna2aa
import pandas as pd
import seaborn as sns
import pickle as pkl
import matplotlib.colors as mcolors
from scipy import stats
import os.path
from matplotlib.lines import Line2D
import json
import shutil
#from evaluation_functions import *
from functions_ import *
from plotting import *
from Bio import SeqIO
import matplotlib.patches as patches
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from characterization_from_blast_alignments import *


data_dir = "data/fastq/P03_RL8_AraCLOV2" ## within data_dir, there should be two directories: 1) /references (containing the reference sequence) and 2) /blast/alignments (containing the blast output files)
#P01_DP6_LOV2/" #P02_RL8_LOV2
with open(f"{data_dir}/config.json", "r") as file:
    config = json.load(file)

read_directions = [ "R1"] #read directions that should be considered for the analysis ["R1", "R2"] or ["R1"] or ["R2"]
datatypes = [ "DNA", "AA", "Codons"] # data types that should be considered for the analysis ["DNA", "AA", "Codons"]

roi_startseq = "ttagccacaa".upper() ## LOV2 start # set region of interest, that has to be included in the reads to be considered for the analysis, e.g. LOV2 start site
roi_endseq = "cggccaaa".upper() ## LOV2 end
filter_for_reads_with_roi = True
cut_to_roi = False # if True, the reads will be filtered for the region of interest, if False, the whole read will be considered for the analysis

variant = config["variant"] 
used_Barcodes = config["used_Barcodes"]
Sections = config["Sections"] 
full_amplicon = config["amplicon"]#[2:]
full_amplicon_AA = translate_dna2aa(full_amplicon)
min_coverage = 100
data_type = "AA"



print("############# calculation for", data_type, "#############")
full_reference = full_amplicon if data_type != "AA" else full_amplicon_AA


key_of_interest = "combined" if len(read_directions) > 1 else read_directions[0]

FigFolder = f"{os.getcwd()}/output/{variant}/blast/{key_of_interest}/plots/{data_type}"
if not os.path.exists(FigFolder):
    os.makedirs(FigFolder)


OutputFolder = f"{os.getcwd()}/output/{variant}/blast/{key_of_interest}/enrichments/{data_type}"
if not os.path.exists(OutputFolder):
    os.makedirs(OutputFolder)

############# 

genetic_code = {
'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }

codons = list(genetic_code.keys())


ecoli_pref = { ### codons used for retron library (RL8) construction
            "A": 'GCG',
            "R": 'CGT',
            "N": 'AAC',
            "D": 'GAT',
            "C": 'TGC',
            "Q": 'CAG',
            "E": 'GAA',
            "G": 'GGC',
            "H": 'CAT',
            "I": 'ATT',
            "L": "CTG",
            "K": 'AAA',
            "M": 'ATG',
            "F": "TTT",
            "P": 'CCG',
            "S": 'AGC',
            "T": 'ACC',
            "W": 'TGG',
            "Y": "TAT",
            "V": 'GTG',
}
for Bc in used_Barcodes:
    print("##############", Bc, "##############")
    for Section in Sections: 
        print("##############", Section, "##############")
        blast_stemFilename = f"{variant}_{Bc}_{Section}_Nt_filt_" ## update accordingly, total name should be e.g. f"RL8_BC1_S1_Nt_filt_R1.out", same stem should be used for the reference file and the blast output file

        if not os.path.exists(f"{data_dir}/references/{blast_stemFilename}ref.fasta"):
            print(f"Reference file {data_dir}/references/{blast_stemFilename}ref.fasta does not exist, please check the reference file path")
            exit()