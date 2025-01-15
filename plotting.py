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
                  samplename = "", 
                  FigFolder = None): 
    """
    plot read coverage 

    coverage_df: df with the coverage of each position, e.g by calling variants_df.sum()
    samplename: name of the variant
    FigFolder: folder to save the figure
    """

    plt.bar(coverage_df.columns, coverage_df)
    plt.xlabel("Position")
    plt.ylabel("Read counts")
    plt.title(f'{samplename} read depth')
    if FigFolder:
        if not os.path.exists(FigFolder):
            os.makedirs(FigFolder)
        plt.savefig(f'{FigFolder}/{samplename}_coverage.pdf')
    plt.show()


def plot_mutation_spectrum(data, 
                           samplename="" , 
                           FigFolder = None, 
                           colormap = "viridis"):
    """
    plot mutation spectrum (%) as heatmap

    data = dataframe with the mutagenic spectrum (rows = reference nt, columns = mutated nt), can be calculated using mut_spectrum()
    savepath = folder path to save the figure
    samplename = name of the sample 
    """
    f, ax = plt.subplots(figsize=(6, 6))
    sns.heatmap(data, annot=True, linewidths=.5, ax=ax, vmin = 0, cbar = False, square = True, linecolor = "black", cmap = colormap)
    plt.xlabel('Mutated base (%)', fontsize = 10)
    plt.ylabel('Reference base (%)', fontsize = 10)
    for _, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(.5)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=1, fontsize=10)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=1, fontsize=10)
    plt.title(f"{samplename} mutagenic spectrum", fontsize = 12)

    if FigFolder:
        if not os.path.exists(FigFolder):
            os.makedirs(FigFolder)
        plt.savefig(f"{FigFolder}/{samplename}_mutagenic_spectrum_perc.pdf")
        
    plt.show()
    plt.clf()


def plot_mut_rate_per_pos(data, 
                          ref_seq,
                          data_type = "DNA",
                          samplename = "", 
                          FigFolder = None):
    """
    plot the mutation rate per position 

    data: df with the relative frequency of each Nt/Codon/AA at each position
    ref_seq: reference DNA sequence
    data_type: set to "DNA", "AA" or "Codons"
    """
    data = data.sum(axis = 0)

    ref_seq = ref_seq if data_type != "AA" else translate_dna2aa(ref_seq)

    if data_type == "Codons":
        x_ticklabels = [ref_seq[i:i+3] for i in range(0, len(ref_seq)//3*3, 3)]
    else: 
        x_ticklabels = [ref for ref in ref_seq]

    plt.figure(figsize=(20,2))
    sns.heatmap(pd.DataFrame(data).T, cmap = "viridis", cbar = True, cbar_kws = {"pad": 0.02, "label": "Mutation rate" },linecolor="black", xticklabels=x_ticklabels, yticklabels=False)
    plt.xlabel("Position")
    plt.xticks(rotation = 2,fontsize=6)
    plt.title(samplename)
    if FigFolder:
        if not os.path.exists(FigFolder):
            os.makedirs(FigFolder)
        plt.savefig(f"{FigFolder}/{samplename}_mutation_rate_per_Nt_position.pdf", bbox_inches="tight")
    plt.show()
    plt.clf()


def plot_mutation_enrichment(variants_df, 
                             ref_seq, 
                             samplename = "",
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
    samplename: plot title
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
    plt.title(samplename, fontsize=20)
    plt.xlabel("Sequence", fontsize = 20)

    if FigFolder:    
        if not os.path.exists(FigFolder):
            os.makedirs(FigFolder)
        plt.savefig(f"{FigFolder}/{samplename}_{data_type}_mutation_enrichment.pdf", bbox_inches="tight")

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

        _ , variant_relative_freq = mask_ref_in_variants_df(ref_seq=ref, variant_df=seq_variants, data_type=data_type)

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
        if not os.path.exists(FigFolder):
            os.makedirs(FigFolder)
        plt.savefig(f"{FigFolder}/{Section}_{data_type}_mutation_enrichment_comparison.pdf", bbox_inches="tight")

    plt.show()
    plt.clf()



def compare_mut_enrichement_for_all(read_dict, 
                                    ref_gene, 
                                    Primer_out_of_triplets, 
                                    Barcodes,
                                    Primer_seq, 
                                    Sections = ["S1", "S2", "S3", "S4"], 
                                    xlim_plot = None,
                                    FigFolder = None, 
                                    data_type = "DNA", 
                                    combine_mut_rates = False,
                                    vmin = 0, 
                                    vmax = None, 
                                    samples = ["Mutagenesis_BC1", "NegPosSelection_BC1", "NegPosSelection_BC2", "Mutagenesis_BC2", "NegPosSelection_BC3", "NegPosSelection_BC4"], 
                                    plt_titles = ["Mutagenesis 1", "Neg Selection 1", "Pos Selection 1", "Mutagenesis 3", "Neg Selection  3", "Pos Selection 3"],
                                    plot_coverage = True, 
                                    color_above_vmax_orange = True, 
                                    show_cbar_for_each = False, 
                                    show_plttitles = True, 
                                    cbar_label = "mutation rate", 
                                    show_only_pos = None, 
                                    AApos_xlabelticks = None, 
                                    fig_size = None, 
                                    bias_per_pos = None, 
                                    return_df = False):
    """
    compare mutation enrichment for all sections as heatmaps with coverages and biases plotted below

    read_dict: dictionary with lists of reads per samples, keys should follow this structure: name1_BCx_Sectionx_R1(/R2), e.g. Mutagenesis_BC1_S1_R1
    ref_gene: reference gene sequence
    Primer_out_of_triplets = dictionary with the number of nucleotides at the beginning of the primer seq before a triplet starts, following the structure {Sectionx_fwd : int, Sectionx_rev : int, Sectiony_fwd : int, ... }
    Primer_seq: dictionary with the fwd and rev primer sequences for each section, following the structure {Sectionx_fwd : seq, Sectionx_rev : seq,  ... }
    Barcodes: dictionary with the barcode sequences, following the structure {BCx_fwd : seq, BCx_rev : seq, BCy_fwd : seq, ... }    
    data_type: set to "DNA", "AA" or "Codons"
    Sections: list of sections to compare (should match the keys of read_dict)
    combine_mut_rates: if True, for each sample, the mutation rates are summed up per position
    samples: list of samples to compare (should match the read_dict keys, and include the respective Barcode, following the structure name1_BCx, e.g Mutagenesis_BC1)
    plot_coverage: if True, coverage is plotted below the mutation enrichment heatmap
    color_above_vmax_orange: if True, values above vmax are colored orange
    show_cbar_for_each: if True, a colorbar is shown for each heatmap, otherwise colorbars are shown below 
    show_plttitles: if True, plt_titles are shown above the heatmaps
    show_only_pos: dictionary with the positions to show for each section
    bias_per_pos: dictionary with the bias per position for each section, that should be plotted below
    return_df: if True, df with the mutation rates is returned (if show_only_pos is set, only the positions of interest are returned, if not, a dict with df per section is returned)
    """
    dataType_handler = {"DNA": gather_nt_variants, "Codons": gather_codon_variants, "AA": gather_AA_variants}
    gather_variants = dataType_handler.get(data_type)
    if not gather_variants: 
        print("Data type not found!")
        exit()

    pltsize = len(samples)+1 if plot_coverage else len(samples)
    pltsize = pltsize + 1 if bias_per_pos else pltsize

    if return_df:
        if show_only_pos: 
            final_df = pd.DataFrame()
            position_labels = []
        else: 
            final_df = {}

    if show_only_pos:
        fig = plt.figure(figsize=(20*len(Sections), 20) if not fig_size else fig_size)
        gs = gridspec.GridSpec(pltsize, len(Sections), height_ratios=[1]*pltsize, width_ratios=[len(show_only_pos[s]) for s in Sections])
        axes = {}

        for i in range(pltsize):
            for j in range(len(Sections)):
                ax = fig.add_subplot(gs[i, j])
                axes[(i, j)] = ax

        fig.subplots_adjust(wspace=0.03)
        
    else:
        fig, axes = plt.subplots(pltsize, len(Sections), figsize=(25*len(Sections), 25) if not fig_size else fig_size)
        fig.subplots_adjust(wspace=0.03)

    for s_idx, Section in enumerate(Sections):

        if return_df: 
            enriched_regions_sec_cylce= pd.DataFrame()

        tripl_st = Primer_out_of_triplets[Section+"_fwd"]
        tripl_end = Primer_out_of_triplets[Section+"_rev"]
        ref_section = find_reference_seq(ref_gene=ref_gene,Section = Section, Primer_seq=Primer_seq, Primer_out_of_triplets=Primer_out_of_triplets)
        ref_gene_section = ref_section

        ref_section = ref_section if data_type != "AA" else translate_dna2aa(ref_section)

        for idx, sample in enumerate(samples):

            Bc = sample[sample.index("BC"):sample.index("BC")+3]
            a_seq = read_dict[sample+ f"_{Section}_R1"]
            b_seq = read_dict[sample + f"_{Section}_R2"]

            seq_variants = gather_variants(a_seq=a_seq, b_seq = b_seq, catch_left=Barcodes[f"{Bc}_fwd"]+Primer_seq[Section + "_fwd"][:tripl_st], catch_right=dna_rev_comp(Barcodes[f"{Bc}_rev"]+Primer_seq[Section+"_rev"][:tripl_end]), ref=ref_section, use_forward_read=True, use_rev_read=True)

            seq_variants = pd.DataFrame.from_dict(seq_variants)

            coverage_df = pd.DataFrame(seq_variants.sum())

            _ , variant_relative_freq = mask_ref_in_variants_df(ref_seq=ref_section, variant_df=seq_variants, data_type=data_type)
            
            xlim_plot = xlim_plot if xlim_plot else variant_relative_freq.shape[1]
            variant_relative_freq = variant_relative_freq.iloc[:,:xlim_plot]

            if combine_mut_rates: 
                variant_relative_freq = pd.DataFrame(variant_relative_freq.sum(axis = 0)).T

            my_cmap = plt.get_cmap('viridis').copy()

            if color_above_vmax_orange:
                my_cmap.set_over('orange')

            if show_only_pos:
                positions = show_only_pos[Section]
                variant_relative_freq = variant_relative_freq.iloc[:,positions]
                ref_section_start = ref_gene[Primer_out_of_triplets["S1_fwd"]:].index(ref_gene_section)//3
                ref_labels = "".join([ref_section[pos] for pos in positions]) if not AApos_xlabelticks else [AApos_xlabelticks[ref_section_start:][pos] for pos in positions]
                coverage_df = coverage_df.iloc[positions,:]

                if return_df and idx == 0:
                    position_labels.extend(ref_labels)
            
            sns.heatmap(variant_relative_freq, annot=False, ax=axes[idx,s_idx], linecolor="black", cmap=my_cmap,  cbar_kws={'label': cbar_label, "pad": 0.02}, vmin=vmin,vmax = vmax, xticklabels=False if idx != len(samples)-1 else True, yticklabels=True if combine_mut_rates == False else False, cbar=show_cbar_for_each )

            for _, spine in axes[idx,s_idx].spines.items():
                spine.set_visible(True)
                spine.set_linewidth(2)
            axes[idx,s_idx].set_yticklabels(axes[idx,s_idx].get_yticklabels(), rotation=1, fontsize=7)
            
            if show_plttitles: 
                axes[idx,s_idx].set_title(plt_titles[idx], fontsize=15)
            
            axes[idx,s_idx].set_facecolor('gray')
            axes[idx,s_idx].grid(False)

            if idx == len(samples)-1:
                axes[idx,s_idx].set_xticklabels(ref_section[:xlim_plot] if not show_only_pos else ref_labels, rotation=1, fontsize=15 if data_type != "DNA" else 7)
            
            if return_df:
                enriched_regions_sec_cylce = pd.concat([enriched_regions_sec_cylce, variant_relative_freq], axis=0)

        if plot_coverage:
            sns.heatmap(coverage_df.T, ax=axes[pltsize-1,s_idx],square=False, cbar_kws={'label': f"coverage pos selection c3", "pad": 0.02}, vmin=0, yticklabels=False, xticklabels=False, vmax=500, cbar=show_cbar_for_each)
        
        if bias_per_pos:
            spec_cmap = sns.light_palette("black", n_colors=30, reverse=False, as_cmap=True)
            vmin_bias = min(([val for values in bias_per_pos.values() for val in values]))
            vmax_bias = max(([val for values in bias_per_pos.values() for val in values]))

            if show_only_pos:
                positions = show_only_pos[Section]
                biases = [bias_per_pos[Section][pos] for pos in positions] ## filter to pos of interest
            else:
                biases = bias_per_pos[Section]

            sns.heatmap(pd.DataFrame(biases).T, ax=axes[pltsize-2, s_idx], cmap=spec_cmap, square=False, cbar_kws={'label': "chance of codon mutation", "pad": 0.02}, yticklabels=False, xticklabels=False, cbar=show_cbar_for_each, vmin=vmin_bias, vmax=vmax_bias)

        if return_df: 
            if show_only_pos: 
                final_df = pd.concat([final_df, enriched_regions_sec_cylce], axis=1)
            else:
                enriched_regions_sec_cylce.index = plt_titles
                final_df[Section] = enriched_regions_sec_cylce
        
    if not show_cbar_for_each:
        if not vmax: 
            print("Please set vmax to show one colorbar for all")
            exit()
        ## add at the bottom of the figure horizontally a cbar for the relative counts
        cbar_ax = fig.add_axes([0.13, 0.05, 0.15, 0.02])
        cbar = fig.colorbar(axes[0,0].collections[0], cax=cbar_ax, orientation = "horizontal")
        cbar.set_label(cbar_label, fontsize = 25)
        cbar.ax.tick_params(labelsize=20)

        ## add cbar also for coverage and biases
        if plot_coverage:
            cbar_ax = fig.add_axes([0.29, 0.05, 0.15, 0.02])
            cbar = fig.colorbar(axes[pltsize-1,0].collections[0], cax=cbar_ax, orientation = "horizontal")
            cbar.set_label('read depth', fontsize = 25)
            cbar.ax.tick_params(labelsize=20)

        if bias_per_pos: 
            cbar_ax = fig.add_axes([0.45, 0.05, 0.15, 0.02])
            cbar = fig.colorbar(axes[pltsize-2,0].collections[0], cax=cbar_ax, orientation = "horizontal")
            cbar.set_label('chance of codon mutation', fontsize = 25)
            cbar.ax.tick_params(labelsize=20)
    
    if FigFolder:
        name = f"PACE_allSections_{data_type}mutation_enrichment_comparison.pdf" if not show_only_pos else f"PACE_allSections_{data_type}mutation_enrichment_comparison_highMutPos.pdf"

        if not os.path.exists(FigFolder):
            os.makedirs(FigFolder)

        plt.savefig(f"{FigFolder}/{name}", bbox_inches="tight")

    plt.show()
    plt.clf()

    if return_df:
        if show_only_pos:
            final_df.columns = position_labels 
            final_df.index = plt_titles
        return final_df


