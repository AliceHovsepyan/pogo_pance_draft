import os
import json

from Bio.SeqIO import QualityIO
import numpy as np
# import re
from DMS_utils import dna_rev_comp, translate_dna2aa
import pandas as pd
from functions_ import *
from plotting import *
from Bio import SeqIO
from Bio.Seq import Seq
from characterization_from_blast_alignments import *

########## define variables for the analysis, please update accordingly
########## run filter_demultiplex_blast.py before running this script

data_dir = "data/fastq/P02_RL8_LOV2" ## within data_dir, there should be two directories: 1) /references (containing the reference sequence) and 2) /blast/alignments (containing the blast output files)
#P01_DP6_LOV2/" #P02_RL8_LOV2
with open(f"{data_dir}/config.json", "r") as file:
    config = json.load(file)

read_directions = ["R1"] #read directions that should be considered for the analysis ["R1", "R2"] or ["R1"] or ["R2"]
data_type = "AA"
roi_startseq = "ttagccacaa".upper() ## LOV2 start # set region of interest, that has to be included in the reads to be considered for the analysis, e.g. LOV2 start site
roi_endseq = "cggccaaa".upper() ## LOV2 end

filter_for_reads_with_roi = True
cut_to_roi = True # if True, the reads will be filtered for the region of interest, if False, the whole read will be considered for the analysis

variant = config["variant"] 
used_Barcodes = config["used_Barcodes"]
Sections = config["Sections"] 
full_amplicon = config["amplicon"]#[2:]
full_amplicon_AA = translate_dna2aa(full_amplicon)


FigFolder = f"{os.getcwd()}/output/{variant}/blast/plots/"
if not os.path.exists(FigFolder):
    os.makedirs(FigFolder)


OutputFolder = f"{os.getcwd()}/output/{variant}/blast/enrichments/"
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


##### analysis of mutation enrichments

for Bc in used_Barcodes:
    print("##############", Bc, "##############")
    for Section in Sections: 
        print("##############", Section, "##############")
        blast_stemFilename = f"{variant}_{Bc}_{Section}_Nt_filt_" ## update accordingly, total name should be e.g. f"RL8_BC1_S1_Nt_filt_R1.out", same stem should be used for the reference file and the blast output file

        if not os.path.exists(f"{data_dir}/references/{blast_stemFilename}ref.fasta"):
            print(f"Reference file {data_dir}/references/{blast_stemFilename}ref.fasta does not exist, please check the reference file path")
            exit()

        amplicon_seq = str(SeqIO.read(f"{data_dir}/references/{blast_stemFilename}ref.fasta", "fasta").seq)
        amplicon_AA = translate_dna2aa(amplicon_seq)
        reference = amplicon_seq if data_type !="AA" else amplicon_AA

        roi_startidx = amplicon_seq.index(roi_startseq)
        roi_endidx = amplicon_seq.index(roi_endseq) + len(roi_endseq)

        #### calculate mutation enrichment
            
        all_enrichments = {read_dir:{} for read_dir in read_directions}

        for read_dir in read_directions:

            if not os.path.exists(f"{data_dir}/blast/alignments/{blast_stemFilename}{read_dir}.out"):
                print("Blast output file does not exist, please check the blast output file path")
                exit
        
            # Open the blast output file and load it as a dictionary
            print("################",  read_dir,   "################")

            with open(f"{data_dir}/blast/alignments/{blast_stemFilename}{read_dir}.out", "r") as file:
                blast_output = json.load(file)

            blast_alignments = blast_output["BlastOutput2"][0]["report"]["results"]["search"]["hits"].copy()

            ### filter blast alignments for regions that include the region of interest (e.g. LOV2 insertion site) (span at least 10 nucleotides before and after the region)
            print(len(blast_alignments), "total alignments")

            if filter_for_reads_with_roi:
                print("Filtering for reads with region of interest")
                filter_for_region = roi_startidx if read_dir=="R1" else roi_endidx

                blast_alignments = [alignment for alignment in blast_alignments if alignment["hsps"][0]["query_from"] <= filter_for_region-10 and alignment["hsps"][0]["query_to"] >= filter_for_region+10]

                print(len(blast_alignments), "alignments after filtering filtering for reads with region of interest")

            alignments = restructure_alignments(blast_alignments, query_seq=amplicon_seq, read_dir=read_dir)


            #### calculate enrichments
            all_variants, indels_freq,  enrichment_counts, enrichment_relative = characterize_DMS_blast_alignment(alignments, amplicon_seq, data_type=data_type,read_dir=read_dir, exclude_not_covered_regions=False if len(read_directions) > 1 else True)
            all_variants = pd.DataFrame.from_dict(all_variants)

            # cut to region of interest if cut_to_roi is True

            cut_start_idx = roi_startidx if cut_to_roi else 0
            cut_start_idx = cut_start_idx//3 if data_type != "DNA" else cut_start_idx
            cut_end_idx = roi_endidx if cut_to_roi else len(amplicon_seq)
            cut_end_idx = cut_end_idx//3 if data_type != "DNA" else cut_end_idx

            ### add correct column indexes, based on "full_amplicon", if both read directions are considered
            if len(read_directions) > 1: 
                idxs_in_full_amplicon = [idx for idx in range(full_amplicon.index(amplicon_seq), full_amplicon.index(amplicon_seq)+len(amplicon_seq))] if data_type == "DNA" else [idx for idx in range(full_amplicon_AA.index(amplicon_AA), full_amplicon_AA.index(amplicon_AA)+len(amplicon_AA))]
                all_variants.columns = idxs_in_full_amplicon
                indels_freq.columns = [idx for idx in range(full_amplicon.index(amplicon_seq), full_amplicon.index(amplicon_seq)+len(amplicon_seq))] ## always DNA level
                enrichment_counts.columns = idxs_in_full_amplicon
                enrichment_relative.columns = idxs_in_full_amplicon

            ### store enrichments
            all_enrichments[read_dir]["all_variants"] = all_variants.iloc[:,cut_start_idx:cut_end_idx]
            all_enrichments[read_dir]["indels_freq"] = indels_freq.iloc[:,cut_start_idx:cut_end_idx]
            all_enrichments[read_dir]["enrichment_counts"] = enrichment_counts.iloc[:,cut_start_idx:cut_end_idx]
            all_enrichments[read_dir]["enrichment_relative"] = enrichment_relative.iloc[:,cut_start_idx:cut_end_idx]
            print(all_enrichments[read_dir]["all_variants"])

        if len(read_directions) > 1: 

            all_enrichments["combined"] = {}

            ## total variants of R1 and R2
            all_enrichments["combined"]["all_variants"] =  all_enrichments["R1"]["all_variants"] + all_enrichments["R2"]["all_variants"]

            ## total enrichments of R1 and R2
            all_enrichments["combined"]["enrichment_counts"] = all_enrichments["R1"]["enrichment_counts"] + all_enrichments["R2"]["enrichment_counts"]
            all_enrichments["combined"]["enrichment_relative"] = all_enrichments["combined"]["enrichment_counts"] / all_enrichments["combined"]["all_variants"].sum()

            ### combine indels of R1 and R2
            all_enrichments["combined"]["indels_freq"] =  all_enrichments["R1"]["indels_freq"] + all_enrichments["R2"]["indels_freq"]

        ### set reference sequence to the region of interest
        if len(read_directions) > 1:
            reference = reference[cut_start_idx:cut_end_idx] if data_type != "Codons" else reference[cut_start_idx*3:cut_end_idx*3]

        else: ###################### das indexen ist falsch irgendwie 
            if data_type == "Codons": 
                print(cut_start_idx*3, all_enrichments["R1"]["enrichment_relative"].shape[1]*3)
                reference[cut_start_idx*3:all_enrichments["R1"]["enrichment_relative"].shape[1]*3] if read_directions[0] == "R1" else reference[-all_enrichments["R1"]["enrichment_relative"].shape[1]*3:cut_end_idx*3]
                print(reference)
            else:
                reference = reference[cut_start_idx:all_enrichments["R1"]["enrichment_relative"].shape[1]] if read_directions[0] == "R1" else reference[-all_enrichments["R1"]["enrichment_relative"].shape[1]:cut_end_idx]

        ### which key to use for the analysis
        key_of_interest = "combined" if len(read_directions) > 1 else read_directions[0]
        filename = f"{variant}_{Bc}_{Section}_{key_of_interest}_roi{cut_to_roi}_{data_type}"

        ### calculate mutation rates

        average_coverage = all_enrichments[key_of_interest]["all_variants"].sum().sum()/(all_enrichments[key_of_interest]["all_variants"].max().max()*all_enrichments[key_of_interest]["all_variants"].shape[1])*100
        print(f'The illumina paired reads cover on average {average_coverage.round(1)} % of the LOV sequence')

        mut_rate = all_enrichments[key_of_interest]["enrichment_relative"].sum().sum()

        print(f'The mutation rate is estimated to be {mut_rate.round(3)} {data_type} mutations per sequence')
        print(f'If we correct for the coverage, we expect a mutation rate of {round(mut_rate/average_coverage*100,3)} {data_type} mutations per sequence')

        mut_rates_dict = {"calculated_on": data_type, 
                        "coverage": average_coverage,
                        "mut_per_sequence": mut_rate,
                        "mut_per_sequence_coverage_corrected": mut_rate/average_coverage*100
                        }
        
        with open(f'{OutputFolder}/{filename}_mutation_rates.json', 'w') as file:
            file.write(json.dumps(mut_rates_dict, indent=4))

        plot_mutation_enrichment(all_enrichments[key_of_interest]["enrichment_relative"] , ref_seq=reference, samplename=filename, data_type=data_type, FigFolder=FigFolder, vmax=None)
        plt.close()

        ### plot coverages
        coverage_plot( all_enrichments[key_of_interest]["all_variants"].sum(), FigFolder=FigFolder, samplename = filename)
        plt.close()

        ### plot indel frequencies
        if data_type == "DNA":
            fig, axes = plt.subplots(1, figsize=(15,5))
            plt.plot( all_enrichments[key_of_interest]["indels_freq"].columns, all_enrichments[key_of_interest]["indels_freq"].loc["insertion",:], label = "insertion", alpha = 0.5)
            plt.plot( all_enrichments[key_of_interest]["indels_freq"].columns, all_enrichments[key_of_interest]["indels_freq"].loc["deletion",:], label = "deletion", alpha = 0.5)
            if roi_startseq != "" and roi_endseq != "":
                plt.axvline(x=roi_endidx, color='grey', linestyle='--', label = "insertion site")
                plt.axvline(x=roi_startidx, color='grey', linestyle='--')
            plt.legend(frameon = False)
            plt.xlabel("Position")
            plt.ylabel("Frequency")
            plt.title(f"Indel frequency {variant}_{Bc}_{Section}")
            plt.savefig(f"{FigFolder}/{filename}_indel_freq.pdf", bbox_inches='tight')
            plt.savefig(f"{FigFolder}/{filename}_indel_freq.png", bbox_inches='tight')
            plt.close()
        
        ## calculate mutation enrichment
        mut_spec, mut_spec_perc = calc_mut_spectrum_from_enrichment(all_variants, ref_seq=reference, data_type=data_type)

        plot_mutation_spectrum(mut_spec_perc, FigFolder=FigFolder, samplename = filename, data_type=data_type)

        ## save enrichments as csv
        mut_spec_perc.to_csv(f"{OutputFolder}/{filename}_mut_spec.csv")
        all_enrichments[key_of_interest]["enrichment_relative"].to_csv(f"{OutputFolder}/{filename}_enrichment_relative.csv")
        all_enrichments[key_of_interest]["all_variants"].to_csv(f"{OutputFolder}/{filename}_all_variants.csv")
        all_enrichments[key_of_interest]["indels_freq"].to_csv(f"{OutputFolder}/{filename}_indels.csv")
        all_enrichments[key_of_interest]["enrichment_counts"].to_csv(f"{OutputFolder}/{filename}_enrichment_counts.csv")