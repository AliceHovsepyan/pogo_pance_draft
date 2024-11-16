import os
from Bio.SeqIO import QualityIO
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import gzip
import glob
import re
from DMS_utils import dna_rev_comp, translate_dna2aa
import pysam
import pandas as pd
import seaborn as sns
import pickle as pkl
import matplotlib.colors as mcolors
from scipy import stats
import os.path
from matplotlib.lines import Line2D
import json
import shutil



def find(string, value_list):
    indexes = [string.find(letter) for letter in value_list]
    try: 
        ind = min([index for index in indexes if index != -1])
    except:
        ind = 250
    return ind

def read_sequences(variant, catch_left, catch_right, base_dir = os.getcwd(), arbitrary_cutoff_a = False, arbitrary_cutoff_b = False, quality_score = ['!', '"', '#', '$', '%', '&', "'", '(', ')', '*','+', ',', '-', '.', '/', '0', '1', '2', '3', '4', '5'], return_qualities_ids = False):
    """
    read sequences from fastq files while filtering for quality score (read is aborted at first nt with higher error rate than 1%)
    arbitrary_cutoff_a: at which position to arbitrary cut off the forward reads that already went through the quality score filter (= max length of the reads)
    arbitrary_cutoff_b: at which position to arbitrary cut off the backward reads that already went through the quality score filter (= max length of the reads)
    returns list of sequences
    """
    a_sequences = []
    b_sequences = []
    a_qualities = []
    b_qualities = []
    a_ids = []
    b_ids = []

    with open(f'{base_dir}/data/fastq/{variant}_R1_001.fastq', "rt") as a_file, open(f'{base_dir}/data/fastq/{variant}_R2_001.fastq', "rt") as b_file:

        a_reader = QualityIO.FastqGeneralIterator(a_file)
        b_reader = QualityIO.FastqGeneralIterator(b_file)
        
        for total_read, (a, b) in enumerate(zip(a_reader, b_reader)):
                a_id, a_seq, a_qual = a
                b_id, b_seq, b_qual = b
                cutoff_a = find(a_qual, quality_score)
                cutoff_b = find(b_qual, quality_score)

                if arbitrary_cutoff_a and catch_left in a_seq: # cut off a_seq to an (arbitrary) chosen maximum length (=arbitrary_cutoff_a)
                    if cutoff_a > (a_seq.index(catch_left) + arbitrary_cutoff_a):
                        cutoff_a = a_seq.index(catch_left)  + len(catch_left) + arbitrary_cutoff_a 
                
                if arbitrary_cutoff_b and dna_rev_comp(catch_right) in b_seq: 
                    if cutoff_b > (b_seq.index(dna_rev_comp(catch_right)) + arbitrary_cutoff_b):
                        cutoff_b =b_seq.index(dna_rev_comp(catch_right))+ len(catch_right) + arbitrary_cutoff_b

                a_sequences.append(a_seq[:cutoff_a])
                a_qualities.append(a_qual[:cutoff_a])
                b_sequences.append(b_seq[:cutoff_b])
                b_qualities.append(b_qual[:cutoff_b])
                a_ids.append(a_id)
                b_ids.append(b_id)
        print("total reads", total_read+1)

    if return_qualities_ids:
        return a_sequences, b_sequences, a_qualities, b_qualities, a_ids, b_ids  
    else: 
        return a_sequences, b_sequences



def read_filtering(a_seqs, b_seqs,ref_gene, catch_left , catch_right, n_mut_treshold = 10 ): 
    """
    filter out reads with more than n_mut_treshold mutations to get rid of reads with indels that lead to frameshifts and scew the results (as we are here only interested in mutations)
    """
    print("total forward reads before filtering", sum([a_Seq != "" for a_Seq in a_seqs]))
    print("total reverse reads before filtering", sum([b_Seq != "" for b_Seq in b_seqs]))

    a_sequences = []
    b_sequences = []

    for a_seq, b_seq in zip(a_seqs, b_seqs):
            if catch_left in a_seq:
                index = a_seq.index(catch_left) + len(catch_left)
                gene_a = a_seq[index:]
                total_muts_a = sum([ref_gene[idx] != gene_a[idx] for idx in range(len(gene_a))])
                if total_muts_a <= n_mut_treshold: 
                    a_sequences.append(a_seq)
                else: 
                    a_sequences.append("")
            else: 
                a_sequences.append("")
                     
            if dna_rev_comp(catch_right) in b_seq:
                index = b_seq.index(dna_rev_comp(catch_right)) + len(catch_right)
                gene_b = dna_rev_comp(b_seq[index:(len(b_seq)-index)//3*3+index])
                total_muts_b = sum([ref_gene[::-1][idx] != gene_b[::-1][idx] for idx in range(len(gene_b))])

                if total_muts_b <= n_mut_treshold:
                    b_sequences.append(b_seq)
                else:
                    b_sequences.append("")
            else: 
                b_sequences.append("")

    print("total forward reads after filtering", sum([a_seq != "" for a_seq in a_sequences]))
    print("total reverse reads after filtering", sum([b_seq != "" for b_seq in b_sequences]))
    return a_sequences, b_sequences


def gather_AA_variants(a_seq, b_seq, ref_prot,catch_left, catch_right,  use_backward_read=True, use_forward_read = True ):
    """
    returns a dictionary with the counts of each amino acid at each position
    n_mut_treshold: maximum number of mutations allowed (to filter out reads with indels leading to frameshifts that scew the results, as we are only interested in mutations)

    """
    mutation_dict = {}
    
    for idx in range(len(ref_prot)):
        mutation_dict[idx] = {'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 
                            'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 
                            'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 
                            'W':0, 'Y':0, '*':0, 'wt':0}
    for a_seq, b_seq in zip(a_seq, b_seq):
        if use_forward_read:
            if catch_left in a_seq:
                index = a_seq.index(catch_left) + len(catch_left)
                gene_a = a_seq[index:]
                tr_a = translate_dna2aa(gene_a)

                for idx, pos in enumerate(tr_a):
                    mutation_dict[idx][pos] += 1

        if use_backward_read: 
            if dna_rev_comp(catch_right) in b_seq:
                index = b_seq.index(dna_rev_comp(catch_right)) + len(catch_right)
                gene_b = dna_rev_comp(b_seq[index:(len(b_seq)-index)//3*3+index])
                tr_b = translate_dna2aa(gene_b)
                tr_b = tr_b[::-1]
            
                for idx, pos in enumerate(tr_b):
                    mutation_dict[len(ref_prot)-idx-1][pos] += 1
    return mutation_dict

def gather_codon_variants(a_seq, b_seq, ref_gene ,catch_left, catch_right, codons , use_backward_read= True,use_forward_read = True, ):
    """
    returns a dictionary with the counts of each codon at each position
    """
    mutation_dict = {}
    gene_len = len(ref_gene)
    
    for idx in range(0, gene_len//3):
        mutation_dict[idx] = {codon: 0 for codon in codons}

    for a_seq, b_seq in zip(a_seq, b_seq):
        if use_forward_read: 
            if catch_left in a_seq:
                index = a_seq.index(catch_left) + len(catch_left)
                gene_a = a_seq[index:]

                for i in range(0,len(gene_a)//3*3,3): # triplets
                    mutation_dict[i//3][gene_a[i:i+3]] += 1

        if use_backward_read:
            if dna_rev_comp(catch_right) in b_seq:
                index = b_seq.index(dna_rev_comp(catch_right)) + len(catch_right)
                gene_b = dna_rev_comp(b_seq[index:(len(b_seq)-index)//3*3+index])
             
                codons = [gene_b[i:i+3] for i in range(0,len(gene_b),3)]
                if len(gene_b) >= 3:
                    for idx in range(len(codons)): # triplets 
                        codon = codons[-idx-1]
                        if codon:
                            mutation_dict[len(mutation_dict)-idx-1][codon] += 1 # start from the end and update the mutation_dict for each triplet

    return mutation_dict



def gather_nt_variants(a_seq, b_seq ,ref_seq, catch_left , catch_right , use_backward_read= True,use_forward_read = True):
    """
    returns a dictionary with the counts of each nt at each position
    """
    mutation_dict = {}
    gene_len = len(ref_seq)
    
    for idx in range(gene_len):
        mutation_dict[idx] = {'A':0, 'T':0, 'G':0, 'C':0}

    for a_seq, b_seq in zip(a_seq, b_seq):
        if use_forward_read: 
            if catch_left in a_seq:
                index = a_seq.index(catch_left) + len(catch_left)
                gene_a = a_seq[index:]

                for idx, pos in enumerate(gene_a):
                    mutation_dict[idx][pos] += 1

        if use_backward_read:
            if dna_rev_comp(catch_right) in b_seq:
                index = b_seq.index(dna_rev_comp(catch_right)) + len(catch_right)
                gene_b = dna_rev_comp(b_seq[index:(len(b_seq)-index)//3*3+index])
                gene_b = gene_b[::-1]
                for idx, pos in enumerate(gene_b):
                    mutation_dict[gene_len-idx-1][pos] += 1

    return mutation_dict


def process_reads(AA_sequence,use_backward_read = True, use_forward_read = True, arbitrary_cutoff_a = False, arbitrary_cutoff_b= False, variants = None, filter_for_n_mut = True, n_mut_treshold=10, base_dir = os.getcwd()):
    """
    process reads for given variants
    use_backward_read: whether or not to use the backward read
    arbitrary_cutoff: where to cut off the forward sequence (maximum length of the reads, otherwise the cutoff is determined by the quality score = 1% error rate)
    if variants = None, all variants stored in the fastq folder are processed
    """
    variants_dict = {}
    path = f'{base_dir}/data/fastq'
    filenames = glob.glob(f'{path}/*')

    if variants is not None: # filter filenames for given variants
        filenames = [path for path in filenames if any(variant in path for variant in variants)]

    for name in filenames: 
        if '_R1' in name:
            name = name.split('/')[-1].split('_R')[0]
            f1 = name
            a_seq, b_seq = read_sequences(f1, arbitrary_cutoff_a = arbitrary_cutoff_a, arbitrary_cutoff_b=arbitrary_cutoff_b)
            if filter_for_n_mut:
                a_seq, b_seq = read_filtering(a_seq, b_seq, n_mut_treshold = n_mut_treshold)
            variants_dict[name] = {}
            variants_dict[name] = get_variants(a_seq,b_seq,use_backward_read=use_backward_read,use_forward_read=use_forward_read)

            print(f'Done: {name}')

        # with open(f'{path}/{variant}_variants.pickle', 'wb') as handle:
        #     pkl.dump(variants_dict, handle)
    return variants_dict

def get_variants(a_seq,b_seq, ref_prot, ref_gene ,catch_right , codons, catch_left , use_backward_read=True,use_forward_read=True):
    
    variants_dict = {}
    variants_dict["AA"] = gather_AA_variants(a_seq, b_seq, use_backward_read=use_backward_read, use_forward_read=use_forward_read, catch_right=catch_right, catch_left=catch_left, ref_prot=ref_prot)
    variants_dict["DNA"] = gather_nt_variants(a_seq, b_seq, use_backward_read=use_backward_read, use_forward_read=use_forward_read, catch_right=catch_right, catch_left=catch_left, ref_seq=ref_gene)
    variants_dict["Codons"] = gather_codon_variants(a_seq, b_seq, codons = codons, use_backward_read=use_backward_read, use_forward_read=use_forward_read, catch_right=catch_right, catch_left=catch_left, ref_gene=ref_gene)

    return variants_dict
