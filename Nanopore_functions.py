import pandas as pd 
import pysam
import os
from Bio import SeqIO
from utils import translate_dna2aa
from functions_ import mask_ref_in_variants_df
import numpy as np


def characterize_DMS_Nanopore(aligned_reads, ref, data_type = "AA"):
    """
    Function to characterize the DMS alignments, by counting the number of insertions, deletions and substitutions per position

    args:
    DMS_alignments: dict, with the sequences of insert that is mutated 
    ref: str, reference DNA sequence 
    data_type: str, "AA", "DNA" or "Codons
    read_dir: str, "R1" or "R2"
    cut_to_same_start: bool, if True, the sequences are cut to the same start position, otherwise, start and end positions have to be provided as query_from and query_to keys

    returns:
    all_variants: dict, with the counts of the variants per position
    indels_freq: pd.DataFrame, with the counts of insertions and deletions per position (normalized to # alignments)
    enrichment_counts: pd.DataFrame, with the counts of the variants per position, with the reference sequence masked
    enrichment_relative: pd.DataFrame, with the relative counts of the variants per position, with the reference sequence masked

    """

    all_variants = {}
    seq_with_off_target_indels = 0
    included_seq = 0
    indels = pd.DataFrame(columns = range(len(ref)), index = ["insertion"], data = 0)

    reference = translate_dna2aa(ref) if data_type == "AA" else ref
    if data_type == "Codons": 
        reference = [reference[i:i+3] for i in range(0, len(reference), 3)]

    if data_type == "AA":
        for idx in range(len(reference)):
            all_variants[idx] = {'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 
                                    'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 
                                    'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 
                                    'W':0, 'Y':0, '*':0, 'X':0} ## X: are triplets with missing nucleotides (i.e. in the aligned seq, "-" is present), that would lead to a frameshift
    elif data_type == "DNA": 
        for idx in range(len(reference)):
            all_variants[idx] = {'A':0, 'C':0, 'G':0, 'T':0, "-":0}
    
    else:
        codons = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT', 'X']
        for idx in range(len(reference)): 
            all_variants[idx] = {codon: 0 for codon in codons}


    for alignment in aligned_reads:
        
        if "-" in alignment: ## track indels
            seq_with_off_target_indels += 1
            shift = 0 ## here, we count the shift of the position compared to the reference, that occurs if there is an insertion in the qseq 

            for idx,nt in enumerate(alignment): 
                pos = idx - shift # adjust for the shift in the index, due to prior insertions

                if nt == "-":
                    indels.loc["insertion", pos] += 1
                    shift += 1 ## to correct for the shift in the index, due to the insertion

        if data_type == "Codons":
            alignment = [alignment[i:i+3] for i in range(0, len(alignment)//3*3, 3)]

        elif data_type == "AA":
            alignment = translate_dna2aa(alignment)
        
        for idx, variant in enumerate(alignment): 
            if "-" in variant:
                variant = "X"
            all_variants[idx][variant] += 1
       
    indels_freq = indels/len(aligned_reads)
    print(seq_with_off_target_indels, "sequences have off target indels")
    print(len(aligned_reads), "sequences are included in the enrichment analysis")

    all_variants = pd.DataFrame.from_dict(all_variants)

    enrichment_counts, enrichment_relative = mask_ref_in_variants_df(variant_df=all_variants, ref_seq=reference if data_type=="AA" else ref, data_type=data_type)
    

    return all_variants, enrichment_counts,enrichment_relative, indels_freq

import re

def read_cleaning_(input_folder, ref, cut_n_bases_from_start=48):
    """
    process aligned reads so that they are forced to be in frame, i.e. if there is a deletion in the read, a "-" is added, if there is a insertion in the read, the position is skipped. 
    This is necessary due to the high error/indel rate of Nanopore sequencing 
    Also makes sure that all reads start at the same position (cut_n_bases_from_start)

    args: 
    input_folder: folder with .bam files that should be processed
    ref: reference sequence, to which the reads were aligned
    cut_n_bases_from_start: number of positions, that should be cut from the reference start (also the reads), so that all of the reads start (in frame) at the same position (if they do align to this position, the reads will be thrown out) --> make sure that after cutting these many bases from the ref start, the read is in frame

    returns: 
    1. a list of all processed reads, 2. a dataframe of indels that are present in the reads after alignment 3. a list of the qualities corresponding to the reads (insertions are skipped, deletions are given empty strings)
    """
    bam_files = [f for f in os.listdir(input_folder) if f.endswith('.bam')]
    all_reads = []
    indels = pd.DataFrame(columns=range(len(ref)), index=["I", "D"], data=0)
    all_qualities = []

    cigar_pattern = re.compile(r"(\d+)([MIDNSHP=X])")  # Regex to extract (length, operation)

    for file_nr, bamfile_name in enumerate(bam_files):
        bam_path = os.path.join(input_folder, bamfile_name)
        bamfile = pysam.AlignmentFile(bam_path, "rb")

        print("Status:", file_nr + 1, "/", len(bam_files), "done")

        for read in bamfile.fetch():
            if read.is_unmapped or read.query_sequence is None:
                print(f"Skipping read {read.query_name}")
                continue

            alignment_start = read.reference_start
            seq = read.query_sequence
            qualitities = read.query_qualities
            refined_qualities = []
            refined_seq_list = []
            ref_pos = 0  # Reference position in the read

            cigar_operations = [(int(length), op) for length, op in cigar_pattern.findall(read.cigarstring)]

            query_pos = 0  # Track position in the read sequence
            ref_pos = alignment_start  # Start position in reference

            for length, operation in cigar_operations:
                if operation == "M":  # Matches (or mismatches)
                    refined_seq_list.extend(seq[query_pos:query_pos + length])
                    refined_qualities.extend(qualitities[query_pos:query_pos + length])
                    query_pos += length
                    ref_pos += length
                elif operation == "I":  # Insertion (extra bases in read)
                    indels.loc["I", ref_pos] += 1# length
                    query_pos += length
                elif operation == "D":  # Deletion (missing bases in read)
                    refined_qualities.extend([""] * length)
                    refined_seq_list.extend("-" * length)
                    indels.loc["D",ref_pos] += 1#length
                    ref_pos += length
                elif operation == "N":  # Skipped region in reference
                    ref_pos += length
                elif operation == "S":  # Soft clipping (ignored bases at ends)
                    query_pos += length
                elif operation == "H":  # Hard clipping (ignored bases, not in read)
                    continue
                elif operation == "P":  # Padding (shouldn't appear in nanopore data)
                    continue

            refined_seq = "".join(refined_seq_list)

            # Cut off bases from the start if needed
            if alignment_start < cut_n_bases_from_start:
                cut_start = cut_n_bases_from_start - alignment_start
                refined_seq = refined_seq[cut_start:]
                refined_qualities = refined_qualities[cut_start:]

                all_reads.append(refined_seq)
                all_qualities.append(refined_qualities)

        print(f"Processed {bamfile_name}")

    print("Total reads:", len(all_reads))
    return all_reads, indels, all_qualities



def get_linker_regions(input_folder, ref, cut_site_seq_left, cut_site_seq_right, left_linker_region_len, right_linker_region_len):
    """ 
    Function to extract the linker regions from the Nanopore reads

    args:
    input_folder: str, path to the folder with the Nanopore bam files
    ref: str, reference sequence
    cut_site_seq_left: seq of the left linker (start of the insert)
    cut_site_seq_right: seq of the right linker (end of the insert)
    left_linker_region_len: int, length of the left linker region (i.e. the number of bases to extract before the left linker pos) = len left linker read
    righ_linker_region_len: int, length of the right linker region (i.e. the number of bases to extract after the right linker pos) = len right linker read

    returns:
    all_left_linkers: list, with the left linker regions
    all_right_linkers: list, with the corresponding right linker regions

    """
    bam_files = [f for f in os.listdir(input_folder) if f.endswith('.bam')]
    all_left_linkers = {}
    all_right_linkers = {}
    indels = pd.DataFrame(columns=range(len(ref)), index=["I", "D"], data=0)
    #all_qualities = []
    left_linker_excluded = 0
    right_linker_excluded = 0
    cigar_pattern = re.compile(r"(\d+)([MIDNSHP=X])")  # Regex to extract (length, operation)

    id_nr = 0

    for file_nr, bamfile_name in enumerate(bam_files):
        bam_path = os.path.join(input_folder, bamfile_name)
        bamfile = pysam.AlignmentFile(bam_path, "rb")

        print("Status:", file_nr + 1, "/", len(bam_files), "done")

        for read in bamfile.fetch():
            if read.is_unmapped or read.query_sequence is None:
                print(f"Skipping read {read.query_name}")
                continue

            
            alignment_start = read.reference_start
            seq = read.query_sequence
            qualitities = read.query_qualities
            refined_seq_list = []
            refined_ref_list = []
            ref_pos = 0  # Reference position in the read

            cigar_operations = [(int(length), op) for length, op in cigar_pattern.findall(read.cigarstring)]

            query_pos = 0  # Track position in the read sequence
            ref_pos = alignment_start  # Start position in reference

            for length, operation in cigar_operations:
                if operation == "M":  # Matches (or mismatches)
                    refined_seq_list.extend(seq[query_pos:query_pos + length])
                    #refined_qualities.extend(qualitities[query_pos:query_pos + length])
                    refined_ref_list.extend(ref[ref_pos:ref_pos+length])
                    query_pos += length
                    ref_pos += length

                elif operation == "I":  # Insertion (extra bases in read)
                    indels.loc["I", ref_pos] += 1# length
                    refined_seq_list.extend(seq[query_pos:query_pos + length])
                    refined_ref_list.extend("-"*length)
                    query_pos += length
                elif operation == "D":  # Deletion (missing bases in read)
                    #refined_qualities.extend([""] * length)
                    refined_seq_list.extend("-" * length)
                    indels.loc["D",ref_pos] += 1#length
                    refined_ref_list.extend(ref[ref_pos:ref_pos+length])
                    ref_pos += length
                elif operation == "N":  # Skipped region in reference
                    ref_pos += length
                elif operation == "S":  # Soft clipping (ignored bases at ends)
                    query_pos += length
                elif operation == "H":  # Hard clipping (ignored bases, not in read)
                    continue
                elif operation == "P":  # Padding (shouldn't appear in nanopore data)
                    continue
            
            refined_seq = "".join(refined_seq_list)
            refined_ref = "".join(refined_ref_list)

            cut_site_left = refined_ref.find(cut_site_seq_left) 
            cut_site_right = refined_ref.find(cut_site_seq_right) 

            if cut_site_left != -1 and cut_site_right != 1: ## if -1, there are insertions in start of LOV2, thus seq not in ref seq and we do not include these seq
                all_left_linkers["id"+str(id_nr)] = {"hseq" : refined_seq[cut_site_left-left_linker_region_len:cut_site_left], 
                                                     "qseq" : refined_ref[cut_site_left-left_linker_region_len:cut_site_left]}
                
                cut_site_right = cut_site_right + len(cut_site_seq_right)
                all_right_linkers["id"+str(id_nr)] = {"hseq": refined_seq[cut_site_right:cut_site_right+right_linker_region_len],
                                                 "qseq": refined_ref[cut_site_right:cut_site_right+right_linker_region_len]}
                
            else: 
                left_linker_excluded +=1
                right_linker_excluded +=1

            id_nr += 1
        print(f"Processed {bamfile_name}")

    print("Total left linkers:", sum([l != "" for l in all_left_linkers]))
    print("Total right linkers:", sum([l != "" for l in all_right_linkers]))
    print(left_linker_excluded, "left linkers are excluded")
    print(right_linker_excluded, "right linkers are excluded")
    return all_left_linkers, all_right_linkers




def get_genotype_dict_from_AAseqs(all_Aas, ref_AAseq, ref_aa_annot, not_masked_positions = None, combined = False): 
    """
    Function to get the genotype dictionary from the AA sequences of the reads

    all_Aas: list of AA sequences of the reads
    ref_AAseq: reference AA sequence
    ref_aa_annot: list of the reference AA annotations
    not_masked_positions: list of positions that should be considered
    combined: whether or not the kind of mutation should be considered when calculating the genotypes

    returns: 
    genotypes: dictionary with the genotypes and their counts
    """
    genotypes = {}

    for read in all_Aas:
        variant = []
        for idx, Aa in enumerate(read): 
            if not_masked_positions and idx not in not_masked_positions: ## if the position is masked, skip
                continue
            else:
                if Aa != "X" and Aa != ref_AAseq[idx]: 
                    observed_mut = f"{ref_aa_annot[idx]}{Aa}" if not combined else ref_aa_annot[idx]
                    variant.append(observed_mut)
            
        if len(variant) == 0: 
            variant = "WT"
        else: 
            variant = "_".join(variant)

        genotypes[variant] = genotypes.get(variant, 0) + 1
        
    return genotypes
        