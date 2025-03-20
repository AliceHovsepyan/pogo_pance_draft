import pandas as pd 
import pysam
import os
from Bio import SeqIO
from utils import translate_dna2aa
from functions_ import mask_ref_in_variants_df
import numpy as np



# def read_cleaning(input_folder,  ref, cut_n_bases_from_start = 9):

#     bam_files = [f for f in os.listdir(input_folder) if f.endswith('.bam')]
#     all_reads = []
#     #all_ids = []
#     indels = pd.DataFrame(columns = (range(len(ref))), index = ["I", "D"], data = 0)

#     for file_nr, bamfile_name in enumerate(bam_files):

#         bam_path = os.path.join(input_folder, bamfile_name)  
#         bamfile = pysam.AlignmentFile(bam_path, "rb") 
#         # get_aligned_pairs(self, matches_only=False, with_seq=False, with_cigar=False)
#         # a list of aligned read (query) and reference positions.
#         # Each item in the returned list is a tuple consisting of the 0-based offset from the start of the read sequence followed by the 0-based reference position.
#         # For inserts, deletions, skipping either query or reference position may be None.
#         # For padding in the reference, the reference position will always be None.
#         print("Status:", file_nr, "/", len(bam_files), "done")
#         # Iterate over reads in the SAM file
#         for read in bamfile.fetch():
#             if read.is_unmapped or read.query_sequence is None:
#                 print(f"Skipping read {read.query_name}")
#                 continue 

#             alignment_start = read.reference_start 
#             seq = read.query_sequence
#             refined_seq_list = []

#             for cigar_tuple in read.get_aligned_pairs(with_seq=True, with_cigar=True): 
#                 if cigar_tuple[-1].value == 0: ## MATCH
#                     refined_seq_list.append(seq[cigar_tuple[0]])

#                 elif cigar_tuple[-1].value == 1: ## INSERTION
#                     indels.loc["I", len(refined_seq_list)] += 1

#                 elif cigar_tuple[-1].value == 2: ## DELETION
#                     refined_seq_list.append("-")
#                     indels.loc["D", cigar_tuple[1]] += 1

#                 # elif cigar_tuple[-1].value == 3: ## SKIPPED REGION i.e. "deletion"
#                 #     refined_seq_list.append("-")
#                 #     indels.loc["N", cigar_tuple[1]] += 1
                
#                 # elif cigar_tuple[-1] == 4: ## SOFT CLIPPING
#                 #     continue ## skip the soft clipped bases

#                 # elif cigar_tuple[-1] == 5: ## HARD CLIPPING
#                 #     continue
                    
#                 # elif cigar_tuple[-1] == 6: ## PADDING
#                 #     print("Padding")
                
#                 # elif cigar_tuple[-1] == 7: ## SEQ MATCH
#                 #     refined_seq_list.append( seq[cigar_tuple[0]])

#                 # elif cigar_tuple[-1] == 8: ## SEQ MISMATCH
#                 #     refined_seq_list.append( seq[cigar_tuple[0]])

#                 # elif cigar_tuple[-1] == 9: ## CBACK
#                 #     print("CBACK")
#                 refined_seq = "".join(refined_seq_list)


#             if alignment_start < cut_n_bases_from_start: 
#                 cut_start = cut_n_bases_from_start-alignment_start
#                 refined_seq = refined_seq[cut_start:]
#                 all_reads.append(refined_seq)
#                 #all_ids.append(read.query_name)
                
#         print(f"Processed {bamfile_name}")

#     print("total reads:", len(all_reads))

#     return all_reads, indels



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

    ref = translate_dna2aa(ref) if data_type == "AA" else ref

    if data_type == "AA":
        for idx in range(len(ref)):
            all_variants[idx] = {'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 
                                    'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 
                                    'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 
                                    'W':0, 'Y':0, '*':0, 'X':0} ## X: are triplets with missing nucleotides (i.e. in the aligned seq, "-" is present), that would lead to a frameshift
    elif data_type == "DNA": 
        for idx in range(len(ref)):
            all_variants[idx] = {'A':0, 'C':0, 'G':0, 'T':0, "-":0}

    for alignment in aligned_reads:
        
        if "-" in alignment: ## track indels
            seq_with_off_target_indels += 1
            shift = 0 ## here, we count the shift of the position compared to the reference, that occurs if there is an insertion in the qseq 

            for idx,nt in enumerate(alignment): 
                pos = idx - shift # adjust for the shift in the index, due to prior insertions

                if nt == "-":
                    indels.loc["insertion", pos] += 1
                    shift += 1 ## to correct for the shift in the index, due to the insertion


        if data_type == "AA":
            alignment = translate_dna2aa(alignment)
        
        for idx, variant in enumerate(alignment): 
            all_variants[idx][variant] += 1
       
    indels_freq = indels/len(aligned_reads)
    print(seq_with_off_target_indels, "sequences have off target indels")
    print(len(aligned_reads), "sequences are included in the enrichment analysis")

    all_variants = pd.DataFrame.from_dict(all_variants)

    enrichment_counts, enrichment_relative = mask_ref_in_variants_df(variant_df=all_variants, ref_seq=ref, data_type=data_type)
    

    return all_variants, enrichment_counts,enrichment_relative, indels_freq

import re

def read_cleaning_(input_folder, ref, cut_n_bases_from_start=48):
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

        # left_linker_region = list(range(left_linker_pos-left_linker_region_len, left_linker_pos))
        # right_linker_region = list(range(right_linker_pos, right_linker_pos+righ_linker_region_len))

        for read in bamfile.fetch():
            if read.is_unmapped or read.query_sequence is None:
                print(f"Skipping read {read.query_name}")
                continue

            
            alignment_start = read.reference_start
            seq = read.query_sequence
            qualitities = read.query_qualities
            #refined_qualities = []
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

            # Cut off bases from the start if needed
            # if alignment_start < cut_n_bases_from_start:
            #     cut_start = cut_n_bases_from_start - alignment_start
            #     refined_seq = refined_seq[cut_start:]
            #     refined_qualities = refined_qualities[cut_start:]

            #     all_reads.append(refined_seq)
            #     all_qualities.append(refined_qualities)

            cut_site_left = refined_ref.find(cut_site_seq_left) 
            cut_site_right = refined_ref.find(cut_site_seq_right) 

            if cut_site_left != -1 and cut_site_right != 1: ## if -1, there are insertions in start of LOV2, thus seq not in ref seq and we do not include these seq
                all_left_linkers["id"+str(id_nr)] = {"hseq" : refined_seq[cut_site_left-left_linker_region_len:cut_site_left], 
                                                     "qseq" : refined_ref[cut_site_left-left_linker_region_len:cut_site_left]}
                
                cut_site_right = cut_site_right + len(cut_site_seq_right)
                all_right_linkers["id"+str(id_nr)] = {"hseq": refined_seq[cut_site_right:cut_site_right+right_linker_region_len],
                                                 "qseq": refined_ref[cut_site_right:cut_site_right+right_linker_region_len]}
                
            else: 
                # all_left_linkers["id"+str(id_nr)] = {"hseq" : "", 
                #                                      "qseq" : refined_ref[cut_site_left-left_linker_region_len:cut_site_left]}
                left_linker_excluded +=1
                # all_right_linkers["id"+str(id_nr)] = {"hseq": "",
                #                                     "qseq": refined_ref[cut_site_right:cut_site_right+right_linker_region_len]}
                right_linker_excluded +=1

            id_nr += 1
        print(f"Processed {bamfile_name}")

    print("Total left linkers:", sum([l != "" for l in all_left_linkers]))
    print("Total right linkers:", sum([l != "" for l in all_right_linkers]))
    print(left_linker_excluded, "left linkers are excluded")
    print(right_linker_excluded, "right linkers are excluded")
    return all_left_linkers, all_right_linkers



def get_linker_variants_for_Nanopore(linker_alignments, wt_linker = "SG", read_dir = "R1"):
    """
    Function to characterize linker variants from the blast alignments

    args: 
    linker_alignments: dict, with the blast algined sequences (hseq, qseq) of the linker
    wt_linker: str, the wildtype linker AA sequence 
    read_dir: str, "R1" or "R2"

    returns:
    linkers: dict, with the counts of the linker variants
    linker_list: list, with all the linker sequences
    """

    frameshifts = 0
    linker_counts = {}
    linker_list = []
    for x in linker_alignments.values():
        linker = ""
        qseq = x["qseq"]
        hseq = x["hseq"]

        is_frameshift_read = (qseq.count("-") - hseq.count("-")) %3 != 0
        ##### Exclude frameshift reads 
        if is_frameshift_read:
            # Insertions (shown as "-" in ref) and deletions that sum up to not multiple of three lead to frameshifts -> exclude these reads
            frameshifts += 1
            continue

        ##### WT sequences
        elif qseq == hseq:  # WT linkers with differences in the rest of the sequence are taken into account below
            linker_counts["wt"] = linker_counts.get("wt", 0) + 1
            linker = "wt"

        ##### Reads with deletions 
        elif (qseq.count("-") - hseq.count("-")) < 0:  # Deletions that are multiple of 3, not leading to frameshifts
                correct_by = 3 - (hseq.count("-") % 3) ## we need to correct for the shift in the index, due to the deletion, since we use a **fixed** length for all hseqs, and not the start of the read, thus, based on this alignment and selection of the region, deletions induce frameshifts (although they may not induce frameshifts in real, due to mapped insertions, i.e. qseq.count("-") - hseq.count("-") %3 )== 0) (this is different to analysis of Illumina from blast alignments, where we always use the whole read, i.e. if we have a deletion, we can just skip the bases in the read to keep in frame (see get_linker_variants_from_blast_alignment)
                del_count = (hseq.count("-") - qseq.count("-")) //3 # Number of deletions in AA level
                hseq_filt = re.sub("-", "", hseq)
                if read_dir == "R2": 
                    linker = translate_dna2aa(hseq_filt)[:len(wt_linker)-(del_count)]  # Linker shortened by 3 Nts = 1 AA
                else: 
                    hseq_filt = hseq_filt[correct_by:]
                    linker = translate_dna2aa(hseq_filt)[-len(wt_linker)+(del_count):]  # Linker shortened by 3 Nts = 1 AA

                linker_counts[linker] = linker_counts.get(linker, 0) + 1

        ###### Reads with substitutions 
        elif (qseq.count("-") - hseq.count("-")) == 0:  # Linker was substituted, but no deletions or insertions present
            hseq_filt = re.sub("-", "", hseq)
            correct_by = 3 - (hseq.count("-") %3) 
            del_count = hseq.count("-")
            if read_dir == "R2":
                if hseq_filt[:len(wt_linker) * 3] == qseq[:len(wt_linker) * 3]:
                    linker_counts["wt"] = linker_counts.get("wt", 0) + 1
                else: 
                    linker = translate_dna2aa(hseq_filt)[:len(wt_linker)]
                    linker_counts[linker] = linker_counts.get(linker, 0) + 1
            else: 
                if hseq_filt[-len(wt_linker) * 3:] == qseq[-len(wt_linker) * 3:]:  
                    # WT linker (but differences in the rest (e.g beginning) of the sequence, thus these did not meet the first criterion)
                    linker_counts["wt"] = linker_counts.get("wt", 0) + 1
                else: 
                    hseq_filt = hseq_filt[correct_by:]
                    linker = translate_dna2aa(hseq_filt)[-len(wt_linker):]
                    linker_counts[linker] = linker_counts.get(linker, 0) + 1

            
        ###### Reads with insertions
        elif (qseq.count("-")- hseq.count("-")) > 0:
            insertion_len = (qseq.count("-") - hseq.count("-")) //3 # AA level  - hseq.count("-") 
            hseq_filt = re.sub("-", "", hseq)
            correct_by = 3 - (hseq.count("-") % 3)
            if read_dir == "R2": 
                linker = translate_dna2aa(hseq_filt)[:len(wt_linker) + insertion_len]
                linker_counts[linker] = linker_counts.get(linker, 0) + 1
            else: 
                hseq_filt = hseq_filt[correct_by:]
                linker = translate_dna2aa(hseq_filt)[-len(wt_linker) - insertion_len:]
                linker_counts[linker] = linker_counts.get(linker, 0) + 1
        #### All other reads
        else:
            print("sequence", hseq, "does not meet any criteria")
        
        linker_list.append(linker)
        
        if linker and linker == "D*RKPAV": 
            print("hseq", hseq)
            print("qseq", qseq)
            print("linker", linker)
    print(frameshifts, "reads excluded due to frameshifts")

    return linker_counts, linker_list
