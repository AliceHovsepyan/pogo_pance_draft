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
from functions_ import *


def coverage_plot(coverage_df, 
                  variant_name = "", 
                  FigFolder = None): 
    """
    plot read coverage 

    coverage_df: df with the coverage of each position, e.g by calling variants_df.sum()
    variant_name: name of the variant, for the title of the plot
    FigFolder: folder to save the figure
    """

    plt.bar(coverage_df.columns, coverage_df)
    plt.xlabel("Position")
    plt.ylabel("Read counts")
    plt.title(f'{variant_name} read depth')
    if FigFolder:
        plt.savefig(f'{FigFolder}/{variant_name}_DNA_coverage.pdf')
    plt.show()


def plot_mutation_enrichment(variants_df, 
                             name, 
                             ref_seq, 
                             data_type = "DNA",
                             FigFolder = None, 
                             cmap = "viridis", 
                             cbar_label = "mutation rate", 
                             vmax = None):
    """
    plot mutation enrichment (heatmap, with ref nts/codons/AA in grey)

    variants_df: df with relative counts of each Nt/AA/Codon (rows) at each position (columns), (call mask_ref_in_variants_dict() prior to set ref Nt/AA/Codon to NA, then shown as grey)
    data_type: set to "DNA", "AA" or "Codons" 
    ref_seq: reference sequence (should be DNA sequence if data_type is "DNA" or "Codons", and AA sequence if data_type is "AA")
    name: plot title
    FigFolder: folder to save the figure
    """

    if data_type in ["DNA", "AA"]:
        seq_pos = list(ref_seq)
    if data_type == "Codons":
        seq_pos = [ref_seq[i:i+3] for i in range(0, len(ref_seq)//3*3, 3)]

    plt.figure(figsize=(30,10))
    sns.reset_defaults()
    
    ax = sns.heatmap(data=variants_df, cmap=cmap, cbar_kws={'label': cbar_label, "pad": 0.02}, yticklabels=True, xticklabels=True, center=0 if cmap == "coolwarm" else None, vmax=vmax, vmin=-vmax if (cmap=="coolwarm" and vmax) else None)

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
    plt.title(name, fontsize=20)
    plt.xlabel("Sequence", fontsize = 20)

    if FigFolder:    
        plt.savefig(f"{FigFolder}/{name}_{data_type}_mutation_enrichment.pdf", bbox_inches="tight")

    plt.show()
    plt.clf()


def compare_mut_enrichement(read_dict, 
                            Section, 
                            ref_gene, 
                            Primer_out_of_triplets, 
                            Barcodes,
                            Primer_seq, 
                            use_rev_read = True, 
                            use_forward_read = True, 
                            xlim_plot = None,
                            FigFolder = None, 
                            data_type = "DNA", 
                            combine_mut_rates =False,
                            vmin = 0, 
                            vmax = None, 
                            samples = ["Mutagenesis_BC1", "NegPosSelection_BC1", "NegPosSelection_BC2", "Mutagenesis_BC2", "NegPosSelection_BC3", "NegPosSelection_BC4"], 
                            plt_titles =["Mutagenesis cycle 1", "Negative selection cycle 1", "Positive selection cycle 1", "Mutagenesis cycle 3", "Negative selection cycle 3", "Positive selection cycle 3"], 
                            plot_coverage = True, 
                            color_above_vmax_orange = True, 
                            cbar_label = "mutation rate", 
                            show_only_pos = None, 
                            fig_size = (20,25)):
    """
    compare mutation enrichment between different mut/selection steps for a given section as heatmap with coverage plotted below

    read_dict: dictionary with reads, following the structure {samplename_Section_R1: reads, samplename_Section_R2: reads, samplename2_S1_R1: reads, ... }
    ref_gene: reference gene sequence
    Primer_seq: dictionary with the fwd and rev primer sequences for each section, following the structure {S1_fwd : seq, S1_rev : seq, S2_fwd : seq, ... }
    Section: Section of interest
    Primer_out_of_triplets = dictionary with the number of nucleotides at the beginning of the primer seq before a triplet starts, following the structure {S1_fwd : int, S1_rev : int, S2_fwd : int, ... }
    Barcodes: dictionary with the barcode sequences, following the structure {BC1_fwd : seq, BC1_rev : seq, BC2_fwd : seq, ... }
    use_forward_read, use_rev_read: whether or not to include the forward read (R1) and/or reverse read (R2) in the analysis (default: True)
    data_type: set to "DNA", "AA" or "Codons"
    combine_mut_rates: if True, for each sample, the mutation rates are summed up per position
    samples: list of samples to compare (should be in the read_dict, and include the respective Barcode, e.g name1_BC1, name2_BC2)
    plot_coverage: if True, coverage is plotted below the mutation enrichment heatmap
    color_above_vmax_orange: if True, values above vmax are colored orange
    show_only_pos: dictionary that contains for each section the positions to show in the heatmap, following the structure {S1: pos, S2: pos, ...}
    """

    dataType_handler = {"DNA": gather_nt_variants, "Codons": gather_codon_variants, "AA": gather_AA_variants}
    gather_variants = dataType_handler.get(data_type)
    if not gather_variants: 
        print("Data type not found!")
        exit()

    tripl_st = Primer_out_of_triplets[Section+"_fwd"]
    tripl_end = Primer_out_of_triplets[Section+"_rev"]
    ref_section = find_reference_seq(ref_gene=ref_gene, Primer_seq=Primer_seq, Section=Section,Primer_out_of_triplets=Primer_out_of_triplets) 

    ref = ref_section if data_type != "AA" else translate_dna2aa(ref_section)

    pltsize = len(samples)+1 if plot_coverage else len(samples)
    fig, axes = plt.subplots(pltsize, 1, figsize=fig_size)
    fig.subplots_adjust(wspace=0.01)

    for idx, sample in enumerate(samples):

        Bc = sample[sample.index("BC"):sample.index("BC")+3]
        a_seq = read_dict[sample + f"_{Section}_R1"]
        b_seq = read_dict[sample+ f"_{Section}_R2"]

        seq_variants = gather_variants(a_seq=a_seq, b_seq = b_seq, catch_left=Barcodes[f"{Bc}_fwd"]+Primer_seq[Section + "_fwd"][:tripl_st], catch_right=dna_rev_comp(Barcodes[f"{Bc}_rev"]+Primer_seq[Section+"_rev"][:tripl_end]), ref=ref, use_forward_read=use_forward_read, use_rev_read=use_rev_read)

        seq_variants = pd.DataFrame.from_dict(seq_variants)
        coverage_df = pd.DataFrame(seq_variants.sum())

        _ , variant_relative_freq = mask_ref_in_variants_dict(ref_seq=ref, variant_df=seq_variants, data_type=data_type)

        if show_only_pos:
            positions = show_only_pos[Section]
            variant_relative_freq = variant_relative_freq.iloc[:,positions]
            ref = "".join([ref[pos] for pos in positions])
            coverage_df = coverage_df.iloc[positions,:]

        xlim_plot = xlim_plot if xlim_plot else variant_relative_freq.shape[1]
        variant_relative_freq = variant_relative_freq.iloc[:,:xlim_plot]

        if combine_mut_rates: 
            variant_relative_freq = pd.DataFrame(variant_relative_freq.sum(axis = 0)).T

        my_cmap = plt.get_cmap('viridis').copy()
        if color_above_vmax_orange:
            my_cmap.set_over('orange')

        sns.heatmap(variant_relative_freq, annot=False, ax=axes[idx], linecolor = "black", cmap = my_cmap,  cbar_kws={'label': cbar_label, "pad": 0.02},vmin=vmin,vmax = vmax,   xticklabels=False if idx != len(samples)-1 else True, yticklabels = True if combine_mut_rates == False else False)

        for _, spine in axes[idx].spines.items():
            spine.set_visible(True)
            spine.set_linewidth(2)
        axes[idx].set_yticklabels( axes[idx].get_yticklabels(), rotation=1, fontsize=5)
        axes[idx].set_title(plt_titles[idx], fontsize = 15)
        axes[idx].set_facecolor('gray')
        axes[idx].grid(False)
        
        if idx == len(samples)-1:
             axes[idx].set_xticklabels(ref[:xlim_plot] , rotation=1, fontsize=7 if data_type != "DNA" else 3)

    if plot_coverage:
        sns.heatmap(coverage_df.T, ax = axes[pltsize-1],square=False, cbar_kws={'label': f"coverage pos selection c3", "pad": 0.02}, vmin = 0, yticklabels= False, xticklabels=False, vmax = 500)
        axes[idx].set_xticklabels(ref[:xlim_plot] , rotation=1, fontsize=7 if data_type != "DNA" else 3)
        
    if FigFolder:
        plt.savefig(f"{FigFolder}/{Section}_{data_type}_mutation_enrichment_comparison.pdf", bbox_inches="tight")

    plt.show()
    plt.clf()