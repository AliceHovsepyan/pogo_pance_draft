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
import matplotlib.gridspec as gridspec




def find(string, 
         value_list
         ):
    
    indexes = [string.find(letter) for letter in value_list]
    try: 
        ind = min([index for index in indexes if index != -1])
    except:
        ind = 250
    return ind


def read_sequences(variant, 
                   catch_left, 
                   catch_right, 
                   base_dir = os.getcwd(), 
                   arbitrary_cutoff_a = False, 
                   arbitrary_cutoff_b = False, 
                   quality_score = ['!', '"', '#', '$', '%', '&', "'", '(', ')', '*','+', ',', '-', '.', '/', '0', '1', '2', '3', '4', '5'], 
                   return_qualities_ids = False):
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


def read_filtering(a_seqs, 
                   b_seqs,
                   ref_gene, 
                   catch_left , 
                   catch_right, 
                   n_mut_treshold = 10 ): 
    """
    filter out reads with more than n_mut_treshold mutations to get rid of reads with indels that lead to frameshifts and scew the results (as we are here only interested in mutations)
    catch_left = left catch sequence
    catch_right = right catch sequence (dna_rev_comp(catch_right) is used to find the sequence in the reverse reads)
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
                # print(ref_gene[::-1])
                # print(gene_b[::-1])
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


def gather_AA_variants(a_seq, 
                       b_seq, 
                       ref_prot,
                       catch_left, 
                       catch_right,  
                       use_rev_read=True, 
                       use_forward_read = True ):
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

        if use_rev_read: 
            if dna_rev_comp(catch_right) in b_seq:
                index = b_seq.index(dna_rev_comp(catch_right)) + len(catch_right)
                gene_b = dna_rev_comp(b_seq[index:(len(b_seq)-index)//3*3+index])
                tr_b = translate_dna2aa(gene_b)
                tr_b = tr_b[::-1]
            
                for idx, pos in enumerate(tr_b):
                    mutation_dict[len(ref_prot)-idx-1][pos] += 1
    return mutation_dict

def gather_codon_variants(a_seq, b_seq, ref_gene ,catch_left, catch_right, codons , use_rev_read= True,use_forward_read = True, ):
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

        if use_rev_read:
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



def gather_nt_variants(a_seq, b_seq ,ref_seq, catch_left , catch_right , use_rev_read= True,use_forward_read = True):
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

        if use_rev_read:
            if dna_rev_comp(catch_right) in b_seq:
                index = b_seq.index(dna_rev_comp(catch_right)) + len(catch_right)
                gene_b = dna_rev_comp(b_seq[index:(len(b_seq)-index)//3*3+index])
                gene_b = gene_b[::-1]
                for idx, pos in enumerate(gene_b):
                    mutation_dict[gene_len-idx-1][pos] += 1

    return mutation_dict


def process_reads(ref_prot, ref_gene,catch_left, catch_right, codons,use_rev_read = True, use_forward_read = True, arbitrary_cutoff_a = False, arbitrary_cutoff_b= False, variants = None, filter_for_n_mut = True, n_mut_treshold=10, base_dir = os.getcwd()):
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
            a_seq, b_seq = read_sequences(f1, arbitrary_cutoff_a = arbitrary_cutoff_a, arbitrary_cutoff_b=arbitrary_cutoff_b, catch_left=catch_left, catch_right=catch_right)
            if filter_for_n_mut:
                a_seq, b_seq = read_filtering(a_seq, b_seq, ref_gene=ref_gene, catch_left= catch_left, catch_right=catch_right,n_mut_treshold = n_mut_treshold)
            variants_dict[name] = {}
            variants_dict[name] = get_variants(a_seq,b_seq,ref_prot = ref_prot, ref_gene = ref_gene, use_rev_read=use_rev_read,use_forward_read=use_forward_read, catch_left=catch_left, catch_right=catch_right, codons = codons)

            print(f'Done: {name}')

        # with open(f'{path}/{variant}_variants.pickle', 'wb') as handle:
        #     pkl.dump(variants_dict, handle)
    return variants_dict

def get_variants(a_seq,b_seq, ref_prot, ref_gene ,catch_right , codons, catch_left , use_rev_read=True,use_forward_read=True):
    
    variants_dict = {}
    variants_dict["AA"] = gather_AA_variants(a_seq, b_seq, use_rev_read=use_rev_read, use_forward_read=use_forward_read, catch_right=catch_right, catch_left=catch_left, ref_prot=ref_prot)
    variants_dict["DNA"] = gather_nt_variants(a_seq, b_seq, use_rev_read=use_rev_read, use_forward_read=use_forward_read, catch_right=catch_right, catch_left=catch_left, ref_seq=ref_gene)
    variants_dict["Codons"] = gather_codon_variants(a_seq, b_seq, codons = codons, use_rev_read=use_rev_read, use_forward_read=use_forward_read, catch_right=catch_right, catch_left=catch_left, ref_gene=ref_gene)

    return variants_dict

def demultiplex_reads(a_seqs, b_seqs,ref_gene, Barcodes , Primer_seq, used_Barcodes, Sections = ["S1", "S2", "S3", "S4"], max_mismatch_primerseq = 1, filter_for_n_mut = True, n_mut_treshold = 20, a_ids = None, b_ids = None):

    read_Dict = {}
    ids_Dict = {}
    ## split the reads into the samples according to Barcode and AraC section -> thereby keeping the forward and reverse reads together
    for Barcode in used_Barcodes: 
        print(Barcode)
        for Section in Sections:
            fwd_BC_Primer_seq = Barcodes[Barcode + "_Fwd"] + Primer_seq[Section+"_fwd_primer"] 
            rev_BC_Primer_seq = Barcodes[Barcode + "_Rev"] + Primer_seq[Section+"_rev_primer"] 
            ### select the reads that contain the forward and reverse BC + primer sequences, thereby allowing for 2 mismatches in the primer sequences but no errors in BCs
            fwd_idxs = [i for i, seq in enumerate(a_seqs) if (seq[:len(Barcodes[Barcode + "_Fwd"])] == Barcodes[Barcode + "_Fwd"] and  sum([sequence!=primer_ref for sequence, primer_ref in zip(seq[len(Barcodes[Barcode + "_Fwd"]):len(fwd_BC_Primer_seq)], Primer_seq[Section+"_fwd_primer"])]) <= max_mismatch_primerseq)]
                        
            rev_idxs = [i for i, seq in enumerate(b_seqs) if (seq[:len(Barcodes[Barcode + "_Rev"] )] == Barcodes[Barcode + "_Rev"] and sum([sequence!=primer_ref for sequence, primer_ref in zip(seq[len(Barcodes[Barcode + "_Rev"]):len(rev_BC_Primer_seq)], Primer_seq[Section+"_rev_primer"])]) <= max_mismatch_primerseq) ]

            indexes = set(fwd_idxs + rev_idxs)
            
            a_seq_Bc_Sec = [a_seqs[i] for i in indexes]
            b_seq_Bc_Sec = [b_seqs[i] for i in indexes]

            if a_ids and b_ids:
                a_ids_Bc_Sec = [a_ids[i].split(" ")[0] for i in indexes]
                b_ids_Bc_Sec = [b_ids[i].split(" ")[0]  for i in indexes]

            ref_seq_Section = ref_gene[ref_gene.index(Primer_seq[Section + "_fwd_primer"]):ref_gene.index(dna_rev_comp(Primer_seq[Section+"_rev_primer"]))+len(Primer_seq[Section+"_rev_primer"])]

            if filter_for_n_mut:
                a_seq_Bc_Sec, b_seq_Bc_Sec = read_filtering(a_seq_Bc_Sec, b_seq_Bc_Sec, catch_left = Barcodes[Barcode + "_Fwd"], catch_right = dna_rev_comp(Barcodes[Barcode + "_Rev"]), n_mut_treshold = n_mut_treshold, ref_gene = ref_seq_Section)

            print(len(a_seq_Bc_Sec), "reads")
            read_Dict[f"{Barcode}_{Section}_R1"] = a_seq_Bc_Sec
            read_Dict[f"{Barcode}_{Section}_R2"] = b_seq_Bc_Sec

            if a_ids and b_ids:
                ids_Dict[f"{Barcode}_{Section}_R1"] = a_ids_Bc_Sec
                ids_Dict[f"{Barcode}_{Section}_R2"] = b_ids_Bc_Sec

    if a_ids and b_ids:
        return  read_Dict, ids_Dict
    else:
        return read_Dict 


def plot_mutation_enrichment(data, name, ref_seq, reverse = False, data_type = "DNA", fig_folder = None, return_df = False, cmap = "viridis", cbar_label = "mutation rate", vmax = None):
    """
    data_type = "DNA", "AA" or "Codons" 
    reference nucleotides/AAs/Codons are shown in grey (set to NA)
    backward: if True, only backward reads are used
    input data should be a dataframe with the relative counts of each nucleotide/AA/Codon at each position
    name = plot title
    """
    #process data
    if data_type == "DNA":
        Nt_order = ['A','C', 'G', 'T']
        data = data.loc[Nt_order]
    elif data_type == "AA": 
        AA_order = ['A','I','L','M','F','W','Y','V','S','T','N','Q','R','H','K','D','E','C','G','P','*']
        data = data.loc[AA_order]
    
    read_len = data.shape[1]
    
    if data_type in ["DNA", "AA"]:
        if reverse: 
            for idx in range(read_len):
                data.loc[ref_seq[::-1][idx], len(ref_seq)-idx-1] = np.nan
            seq_pos = [x for x in ref_seq[-read_len:]] ## for xlabel in plot
            
        else: 
            for idx in range(read_len):
                data.loc[ref_seq[idx], idx] = np.nan
        
            seq_pos = [x for x in ref_seq[:read_len]] ## for xlabel in plot

    elif data_type == "Codons":
        codons = [ref_seq[idx:idx+3] for idx in range(0,len(ref_seq),3)]
        if reverse: 
            for idx in range(read_len):
                data.loc[codons[::-1][idx], len(codons)-idx-1] = np.nan
            seq_pos = [x for x in codons[-read_len:]] ## for xlabel in plot
        
        else: 
            for idx in range(read_len):
                data.loc[codons[idx], idx] = np.nan
            seq_pos = [x for x in codons[:read_len]]## for xlabel in plot

    if return_df: 
        return data
    else: 
        plt.figure(figsize=(30,10))
        sns.reset_defaults()
        
        #sns.set(font_scale =5)
        ax = sns.heatmap(data=data, cmap=cmap, cbar_kws={'label': cbar_label, "pad": 0.02}, yticklabels=True, xticklabels = True, center=0 if cmap == "coolwarm" else None, vmax = vmax, vmin = -vmax if (cmap=="coolwarm" and vmax) else None)
        plt.title(name, fontsize=20)
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(2)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=1, fontsize=10)
        ax.xaxis.set_tick_params(width=2)
        rotation = 90 if data_type == "Codons" else 1
        ax.set_xticklabels(seq_pos, rotation=rotation, fontsize=10)
        ax.yaxis.set_tick_params(width=2)
        ax.set_facecolor('gray')
        ax.grid(False)
        plt.xlabel("sequence", fontsize = 20)

        if fig_folder is not None:    
            plt.savefig(f"{fig_folder}/{name}_mutation_enrichement.pdf", bbox_inches="tight")


        plt.show()
        plt.clf()



def compare_mut_enrichement(read_dict, Section, ref_gene, Primer_out_of_triplets, Barcodes ,Primer_seq , codons, use_rev_read =True, use_forward_read= True, xlim_plot = None,FigFolder = None, data_type = "DNA", combine_mut_rates =False,vmin = 0, vmax =None, variants = ["Mutagenesis_BC1", "NegPosSelection_BC1", "NegPosSelection_BC2", "Mutagenesis_BC2", "NegPosSelection_BC3", "NegPosSelection_BC4"], plt_titles =["Mutagenesis cycle 1", "Negative selection cycle 1", "Positive selection cycle 1", "Mutagenesis cycle 3", "Negative selection cycle 3", "Positive selection cycle 3"], plot_coverage = True, color_above_vmax_red = True, cbar_label = "mutation rate", show_only_pos = None):
    """
    compare mutation enrichment between different mut/selection steps for a given section as heatmap with coverage plotted below
    """

    tripl_st = Primer_out_of_triplets[Section+"_fwd_primer"]
    tripl_end = Primer_out_of_triplets[Section+"_rev_primer"]
    ref_gene_section = ref_gene[ref_gene.index(Primer_seq[Section + "_fwd_primer"][tripl_st:]):ref_gene.index(dna_rev_comp(Primer_seq[Section+"_rev_primer"][tripl_end:]))+len(Primer_seq[Section+"_rev_primer"][tripl_end:])]

    ref_prot_section = translate_dna2aa(ref_gene_section)
    pltsize = len(variants)+1 if plot_coverage else len(variants)

    fig, axes = plt.subplots(pltsize, 1, figsize=(20, 25))
    fig.subplots_adjust(wspace=0.01)

    for idx, variant in enumerate(variants):
        Bc = variant[variant.index("BC"):variant.index("BC")+3]
        a_seq = read_dict[variant + f"_{Section}_R1"]
        b_seq = read_dict[variant+ f"_{Section}_R2"]

        seq_variants = get_variants(a_seq=a_seq, b_seq = b_seq, catch_left=Barcodes[f"{Bc}_Fwd"]+Primer_seq[Section + "_fwd_primer"][:tripl_st],catch_right=dna_rev_comp(Barcodes[f"{Bc}_Rev"]+Primer_seq[Section+"_rev_primer"][:tripl_end]), ref_prot = ref_prot_section, ref_gene = ref_gene_section, use_forward_read=use_forward_read, use_rev_read=use_rev_read, codons = codons)

        seq_variants = pd.DataFrame.from_dict(seq_variants[data_type])
        coverage_df = pd.DataFrame(seq_variants.sum())

        seq_variants = seq_variants/seq_variants.sum()
        ref =  ref_gene_section if data_type!= "AA" else ref_prot_section

        plot_df = plot_mutation_enrichment(data = seq_variants, name = variant, ref_seq = ref , reverse = use_rev_read, data_type = data_type, return_df=True)

        if show_only_pos:
            positions = show_only_pos[Section]
            plot_df = plot_df.iloc[:,positions]
            ref = "".join([ref[pos] for pos in positions])
            coverage_df = coverage_df.iloc[positions,:]


        xlim_plot = xlim_plot if xlim_plot else plot_df.shape[1]
        plot_df = plot_df.iloc[:,:xlim_plot]

        if combine_mut_rates: 
            plot_df = pd.DataFrame(plot_df.sum(axis = 0)).T

        my_cmap = plt.get_cmap('viridis').copy()
        if color_above_vmax_red:
            my_cmap.set_over('orange')

        sns.heatmap(plot_df, annot=False, ax=axes[idx], linecolor = "black", cmap = my_cmap,  cbar_kws={'label': cbar_label, "pad": 0.02},vmin=vmin,vmax = vmax,   xticklabels=False if idx != len(variants)-1 else True, yticklabels = True if combine_mut_rates == False else False)

        for _, spine in axes[idx].spines.items():
            spine.set_visible(True)
            spine.set_linewidth(2)
        axes[idx].set_yticklabels( axes[idx].get_yticklabels(), rotation=1, fontsize=5)
        axes[idx].set_title(plt_titles[idx], fontsize = 15)
        axes[idx].set_facecolor('gray')
        axes[idx].grid(False)
        
        if idx == len(variants)-1:
             axes[idx].set_xticklabels(ref[:xlim_plot] , rotation=1, fontsize=7 if data_type != "DNA" else 3)

    if plot_coverage:
        sns.heatmap(coverage_df.T, ax = axes[pltsize-1],square=False, cbar_kws={'label': f"coverage pos selection c3", "pad": 0.02}, vmin = 0, yticklabels= False, xticklabels=False, vmax = 500)
        axes[idx].set_xticklabels(ref[:xlim_plot] , rotation=1, fontsize=7 if data_type != "DNA" else 3)
        
    if FigFolder:
        plt.savefig(f"{FigFolder}/{Section}_mutation_enrichment_comparison.pdf", bbox_inches="tight")
    plt.show()
    plt.clf()



def compare_mut_enrichement_for_all(read_dict, ref_gene, Primer_out_of_triplets, Barcodes ,Primer_seq , codons, Sections = ["S1", "S2", "S3", "S4"], use_rev_read =True, use_forward_read= True, xlim_plot = None,FigFolder = None, data_type = "DNA", combine_mut_rates =False,vmin = 0, vmax =None, variants = ["Mutagenesis_BC1", "NegPosSelection_BC1", "NegPosSelection_BC2", "Mutagenesis_BC2", "NegPosSelection_BC3", "NegPosSelection_BC4"], plt_titles =["Mutagenesis 1", "Neg Selection 1", "Pos Selection 1", "Mutagenesis 3", "Neg Selection  3", "Pos Selection 3"], plot_coverage = True, color_above_vmax_red = True, show_cbar_for_each=False, show_plttitles = True, cbar_label = "mutation rate", show_only_pos = None, xlabelticks = None):
    """
    compare mutation enrichment for all sections as heatmap with coverage plotted below
    show_only_pos: dictionary with the positions to show for each section
    """

    pltsize = len(variants)+1 if plot_coverage else len(variants)

    if show_only_pos:
        fig = plt.figure(figsize=(20*len(Sections), 20))
        gs = gridspec.GridSpec(pltsize, len(Sections), height_ratios=[1]*pltsize, width_ratios=[len(show_only_pos[s]) for s in Sections])

        axes = {}
        for i in range(pltsize):
            for j in range(len(Sections)):
                ax = fig.add_subplot(gs[i, j])
                axes[(i, j)] = ax
        fig.subplots_adjust(wspace=0.03)
        
    else:
        fig, axes = plt.subplots(pltsize, len(Sections), figsize=(25*len(Sections), 25))
        fig.subplots_adjust(wspace=0.03)

    for s_idx, Section in enumerate(Sections):

        tripl_st = Primer_out_of_triplets[Section+"_fwd_primer"]
        tripl_end = Primer_out_of_triplets[Section+"_rev_primer"]
        ref_gene_section = find_reference_seq(ref_gene=ref_gene,Section = Section, Primer_seq=Primer_seq, Primer_out_of_triplets=Primer_out_of_triplets)

        ref_prot_section = translate_dna2aa(ref_gene_section)

        for idx, variant in enumerate(variants):
            Bc = variant[variant.index("BC"):variant.index("BC")+3]
            a_seq = read_dict[variant+ f"_{Section}_R1"]
            b_seq = read_dict[variant + f"_{Section}_R2"]

            seq_variants = get_variants(a_seq=a_seq, b_seq = b_seq, catch_left=Barcodes[f"{Bc}_Fwd"]+Primer_seq[Section + "_fwd_primer"][:tripl_st],catch_right=dna_rev_comp(Barcodes[f"{Bc}_Rev"]+Primer_seq[Section+"_rev_primer"][:tripl_end]), ref_prot = ref_prot_section, ref_gene = ref_gene_section, use_forward_read=use_forward_read, use_rev_read=use_rev_read, codons = codons)

            seq_variants = pd.DataFrame.from_dict(seq_variants[data_type])

            coverage_df = pd.DataFrame(seq_variants.sum())

            seq_variants = seq_variants/seq_variants.sum()
            ref =  ref_gene_section if data_type!= "AA" else ref_prot_section

            plot_df = plot_mutation_enrichment(data = seq_variants, name = variant, ref_seq = ref , reverse = use_rev_read, data_type = data_type, return_df=True)
            
            xlim_plot = xlim_plot if xlim_plot else plot_df.shape[1]
            plot_df = plot_df.iloc[:,:xlim_plot]

            if combine_mut_rates: 
                plot_df = pd.DataFrame(plot_df.sum(axis = 0)).T

            my_cmap = plt.get_cmap('viridis').copy()
            if color_above_vmax_red:
                my_cmap.set_over('orange')

            if show_only_pos:
                positions = show_only_pos[Section]
                plot_df = plot_df.iloc[:,positions]
                ref_section_start = ref_gene[Primer_out_of_triplets["S1_fwd_primer"]:].index(ref_gene_section)//3
                ref = "".join([ref[pos] for pos in positions]) if not xlabelticks else [xlabelticks[ref_section_start:][pos] for pos in positions]
                coverage_df = coverage_df.iloc[positions,:]
            
            sns.heatmap(plot_df, annot=False, ax=axes[idx,s_idx], linecolor = "black", cmap = my_cmap,  cbar_kws={'label': cbar_label, "pad": 0.02},vmin=vmin,vmax = vmax,   xticklabels=False if idx != len(variants)-1 else True, yticklabels = True if combine_mut_rates == False else False, cbar = show_cbar_for_each )

            for _, spine in axes[idx,s_idx].spines.items():
                spine.set_visible(True)
                spine.set_linewidth(2)
            axes[idx,s_idx].set_yticklabels( axes[idx,s_idx].get_yticklabels(), rotation=1, fontsize=7)
            if show_plttitles: 
                axes[idx,s_idx].set_title(plt_titles[idx], fontsize = 15)
            
            axes[idx,s_idx].set_facecolor('gray')
            axes[idx,s_idx].grid(False)
            # if idx == 0:
            #     axes[idx,s_idx].set_title(Section, fontsize = 30)
            
            if idx == len(variants)-1:
                axes[idx,s_idx].set_xticklabels(ref[:xlim_plot] if not xlabelticks else ref, rotation=1, fontsize=15 if data_type != "DNA" else 7)

        if plot_coverage:
            sns.heatmap(coverage_df.T, ax = axes[pltsize-1,s_idx],square=False, cbar_kws={'label': f"coverage pos selection c3", "pad": 0.02}, vmin = 0, yticklabels= False, xticklabels=False, vmax = 500, cbar = show_cbar_for_each)
            axes[idx,s_idx].set_xticklabels(ref[:xlim_plot] if not xlabelticks else ref, rotation=1, fontsize=15 if data_type != "DNA" else 7)
        
    if not show_cbar_for_each:
        ## add at the bottom of the figure horizontally a cbar for the relative counts
        cbar_ax = fig.add_axes([0.13, 0.05, 0.15, 0.02])
        cbar = fig.colorbar(axes[0,0].collections[0], cax=cbar_ax, orientation = "horizontal")
        cbar.set_label(cbar_label, fontsize = 25)
        cbar.ax.tick_params(labelsize=20)

        # ## ad cbar also for coverage
        cbar_ax = fig.add_axes([0.29, 0.05, 0.15, 0.02])
        cbar = fig.colorbar(axes[pltsize-1,0].collections[0], cax=cbar_ax, orientation = "horizontal")
        cbar.set_label('read depth', fontsize = 25)
        cbar.ax.tick_params(labelsize=20)
        
    if FigFolder:
        name = "PACE_allSections_mutation_enrichment_comparison.pdf" if not show_only_pos else "PACE_allSections_mutation_enrichment_comparison_highMutPos.pdf"
        plt.savefig(f"{FigFolder}/{name}", bbox_inches="tight")
    plt.show()
        

def find_reference_seq(ref_gene, Primer_seq, Section, Primer_out_of_triplets):
    """
    ref_gene = reference gene sequence
    catch_left = left catch sequence
    catch_right = right catch sequence
    Primer_seq = dictionary with primer sequences
    Sections = list of sections
    Primer_out_of_triplets = dictionary with the number of nucleotides at the beginning of the primer seq before a triplet starts
    """ 
    tripl_st = Primer_out_of_triplets[Section+"_fwd_primer"]
    tripl_end = Primer_out_of_triplets[Section+"_rev_primer"]
    ref_gene_section = ref_gene[ref_gene.index(Primer_seq[Section + "_fwd_primer"][tripl_st:]):ref_gene.index(dna_rev_comp(Primer_seq[Section+"_rev_primer"][tripl_end:]))+len(Primer_seq[Section+"_rev_primer"][tripl_end:])]

    return ref_gene_section


## calculate fold change to compare steps/cycles
def calculate_log_FC(read_dictionary, stepA, stepB, Section, BarcodeA, BarcodeB, Primer_out_of_triplets,codons,Primer_seq, Barcodes, ref_gene, data_type="AA", combine_mut_rates = False): 
    """
    all_reads = dictionary with the reads of each step
    stepA, stepB = names of the steps to compare
    normalizeA_by, normalizeB_by = name of the step to normalize the reads by
    """
    ## read all reads 
    stepA_areads = read_dictionary[stepA+"_R1"]
    stepA_breads = read_dictionary[stepA+"_R2"]
    stepB_areads = read_dictionary[stepB+"_R1"]
    stepB_breads = read_dictionary[stepB+"_R2"]

    ref_gene_section = find_reference_seq(ref_gene = ref_gene, Primer_seq = Primer_seq, Section = Section, Primer_out_of_triplets = Primer_out_of_triplets)

    tripl_st = Primer_out_of_triplets[Section+"_fwd_primer"]
    tripl_end = Primer_out_of_triplets[Section+"_rev_primer"]

    ref_prot_section = translate_dna2aa(ref_gene_section)

    ##get all variants 
    stepA_variants = get_variants(a_seq=stepA_areads, b_seq = stepA_breads, catch_left=Barcodes[f"{BarcodeA}_Fwd"]+Primer_seq[Section + "_fwd_primer"][:tripl_st],catch_right=dna_rev_comp(Barcodes[f"{BarcodeA}_Rev"]+Primer_seq[Section+"_rev_primer"][:tripl_end]), ref_prot = ref_prot_section, ref_gene = ref_gene_section, use_forward_read=True, use_rev_read=True, codons = codons)

    stepA_variants = pd.DataFrame.from_dict(stepA_variants[data_type])
    coverage_A = pd.DataFrame(stepA_variants.sum())
    stepA_variants = stepA_variants+1 ## add pseudocount

    ## relative counts
    stepA_variants = stepA_variants/stepA_variants.sum()

    stepB_variants = get_variants(a_seq=stepB_areads, b_seq = stepB_breads, catch_left=Barcodes[f"{BarcodeB}_Fwd"]+Primer_seq[Section + "_fwd_primer"][:tripl_st],catch_right=dna_rev_comp(Barcodes[f"{BarcodeB}_Rev"]+Primer_seq[Section+"_rev_primer"][:tripl_end]), ref_prot = ref_prot_section, ref_gene = ref_gene_section, use_forward_read=True, use_rev_read=True, codons = codons)

    stepB_variants = pd.DataFrame.from_dict(stepB_variants[data_type])
    coverage_B = pd.DataFrame(stepB_variants.sum())
    stepB_variants = stepB_variants+1 ## add pseudocount
    stepB_variants = stepB_variants/stepB_variants.sum()

    ref =  ref_gene_section if data_type!= "AA" else ref_prot_section

    ## calculate Foldchange        

    if combine_mut_rates: 
        for idx in range(len(ref)): ## mask reference AAs/Nts and sum all others at each position --> to then compare the mutation rates
            stepA_variants.loc[ref[idx], idx] = np.nan
            stepB_variants.loc[ref[idx], idx] = np.nan
        stepA_variants = stepA_variants.sum()
        stepB_variants = stepB_variants.sum()

    FC_variants = np.log2(stepA_variants)-np.log2(stepB_variants)    ##pseudocount was already added in the calculation of the relative counts 

    return FC_variants, ref, coverage_A, coverage_B
        
def find_mutated_pos(read_dict, Barcode, Barcodes, Section, ref_gene, Primer_seq, Primer_out_of_triplets, codons, data_type = "AA", cyclename = "Mutagenesis", filter_treshold = 0.05, cov_filter_treshold=0):
    """
    find the positions with a mutation rate above the filter_treshold and the positions with a coverage below the cov_filter_treshold

    read_dict = dictionary with the reads (following this naming convention: {cyclename}_{Barcode}_{Section}_R1:[read1_a, read2_a], {cyclename}_{Barcode}_{Section}_R2: [read1_b, read2_b],...})
    Barcode = name of the barcode
    Barcodes = dictionary with the barcode sequences
    Section = name of the section
    ref_gene = reference gene sequence
    Primer_seq = dictionary with primer sequences
    Primer_out_of_triplets = dictionary with the number of nucleotides at the beginning of the primer seq before a triplet starts
    codons = list of codons
    data_type = "AA" or "DNA"
    cyclename = name of the cycle
    filter_treshold = treshold for the mutation rate
    cov_filter_treshold = treshold for the coverage
    """

    ref_seq_Section = find_reference_seq(ref_gene = ref_gene, Primer_seq = Primer_seq, Section = Section, Primer_out_of_triplets = Primer_out_of_triplets)
    if data_type =="AA":
        ref_prot_Section = translate_dna2aa(ref_seq_Section)

    tripl_st = Primer_out_of_triplets[Section+"_fwd_primer"]
    tripl_end = Primer_out_of_triplets[Section+"_rev_primer"]

    seq_variants = get_variants(read_dict[f"{cyclename}_{Barcode}_{Section}_R1"], read_dict[f"{cyclename}_{Barcode}_{Section}_R2"], catch_left=Barcodes[f"{Barcode}_Fwd"]+Primer_seq[Section + "_fwd_primer"][:tripl_st],catch_right=dna_rev_comp(Barcodes[f"{Barcode}_Rev"]+Primer_seq[Section+"_rev_primer"][:tripl_end]), ref_prot = ref_prot_Section, ref_gene = ref_seq_Section, codons = codons)
    seq_variants = pd.DataFrame.from_dict(seq_variants[data_type])

    coverages = seq_variants.sum()
    ## rates 
    seq_variants = seq_variants/seq_variants.sum(axis = 0)

    ref_seq = ref_seq_Section if data_type == "DNA" else ref_prot_Section
    for idx, ref in enumerate(ref_seq): ## mask reference AAs/Nts and sum all others at each position --> to then compare the mutation rates
        seq_variants.loc[ref, idx] = np.nan

    ## combine mutation rates
    seq_variants = seq_variants.sum(axis = 0)
    low_cov_pos = coverages[coverages<cov_filter_treshold].index 
    high_mut_positions = seq_variants[seq_variants > filter_treshold].index

    return list(high_mut_positions), list(low_cov_pos)