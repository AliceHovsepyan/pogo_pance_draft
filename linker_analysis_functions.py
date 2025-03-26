
import pandas as pd 
import pysam
import os
from Bio import SeqIO
from utils import translate_dna2aa
from functions_ import mask_ref_in_variants_df
import numpy as np
import re

def get_linker_variants(linker_alignments, wt_linker = "SG", read_dir = "R1"):
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
        
        if read_dir == "R1" and len(hseq)%3 != 0: 
            hseq = hseq[(len(hseq)%3):]
            qseq = qseq[(len(qseq)%3):]

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



def rename_right_linkers(linkernames, linkerperc_dict): 
    linkers_perc_filt_corr_names = {}
    linker_renaming = {}

    for startname in linkernames: 
        name = startname
        end_n_del = 0
        start_n_del = 0
        ind = 0
        
        if name.endswith("P"):
            name = name[:-1]
        else: 
            end_n_del += 1
        
        if name.endswith("P") and end_n_del == 0:
            name = name[:-1]
        else: 
            end_n_del += 1
        
        if name.endswith("H") and end_n_del == 0:
            name = name[:-1] 
        else:
            end_n_del += 1

        if name.endswith("L") and end_n_del == 0:
            name = name[:-1]
        else: 
            end_n_del += 1

        if name.startswith("I"): 
            ind +=1
            name = name[1:]
        else: 
            start_n_del += 1
        
        if name.startswith("D"):
            name = name[1:]
            ind +=1
        else:
            start_n_del += 1
        
        if name.startswith("E"):
            name = name[1:]
            ind +=1
        else:
            start_n_del += 1
        
        if name.startswith("A"):
            name = name[1:]
            ind +=1
        else:
            start_n_del += 1
        
        if name.startswith("A"):
            name = name[1:]
            ind +=1
        else:
            start_n_del += 1
        
        if name.startswith("K"):
            name = name[1:]
            ind +=1
        else:
            start_n_del += 1
        
        if end_n_del > 0:
            name = f"{name}(+{end_n_del}del)"
        if start_n_del > 0:
            name = f"(-{start_n_del}del){name}"
        
        print(startname , "->", name)
        linkers_perc_filt_corr_names[name] = linkerperc_dict[startname]
        linker_renaming[startname] = name
        
    return linkers_perc_filt_corr_names, linker_renaming
    

def rename_left_linkers(linkernames, linkerperc_dict): 
    linkers_perc_filt_corr_names = {} 
    linker_renamings = {}
    for startname in linkernames: 
        name = startname
        start_n_del = 0
        end_n_del = 0
        
        if name.endswith("L"): 
            name = name[:-1]
        else:
            end_n_del += 1
        
        if name.startswith("I"):
            name = name[1:]
        else:
            start_n_del += 1
        
        if name.startswith("N"):
            name = name[1:]
        else:
            start_n_del += 1
        
        if name.startswith("E"):
            name = name[1:]
        else:
            start_n_del += 1
        
        if name.startswith("S"):
            name = name[1:]
        else:
            start_n_del += 1
      
        if start_n_del > 0:
            name = f"(-{start_n_del}del){name}"
        if end_n_del > 0:
            name = f"{name}(+{end_n_del}del)"
        
        print(startname , "->", name)
        linkers_perc_filt_corr_names[name] = linkerperc_dict[startname]
        linker_renamings[startname] = name
        
    return linkers_perc_filt_corr_names, linker_renamings
    