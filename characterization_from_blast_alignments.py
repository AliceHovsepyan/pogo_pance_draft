import pandas as pd
import re
from DMS_utils import translate_dna2aa
from functions_ import mask_ref_in_variants_df



def divide_alignments(blast_alignments, cut_site_seq, read_dir="R1"): 
    """ 
    Function to divide the alignments into the linker and insert (LOV2) sequences
    
    args:
    alignments: list of blast alignments 
    cut_site: DNA sequence, of the insert start (if read_dir is R1) or end (if read_dir is R2)
    read_dir: str, "R1" or "R2"

    returns:
    linker_alignments: dict, with the sequences of the linker
    LOV2_alignments: dict, with the sequences of the LOV2 insert
    
    """

    linker_alignments = {}
    LOV2_alignments = {}
    LOV2_start_indel_count = 0


    for alignment in blast_alignments:
        qseq = alignment["hsps"][0]["qseq"].upper()
        hseq = alignment["hsps"][0]["hseq"].upper()
        seq_id = alignment["description"][0]["title"]
        midline = alignment["hsps"][0]["midline"]

        cut_site = qseq.find(cut_site_seq) ## find LOV2 position, we add len(LOV_endseq) if "R2" later, so that we can first filter out reads that not conain the seq of interest (i.e. cut_site = -1) due to insertions at these sites 
        if cut_site != -1: ## if -1, there are insertions in start of LOV2, thus seq not in ref seq and we do not include these seq
            
            if read_dir=="R2":
                cut_site += len(cut_site_seq)

                read_out_of_frame = cut_site % 3  ## correct for out of frame reads, due to different lengths of our reads (only for R2, since R1 is always in frame, begin of the read is cut, but for R2, the end of the read is cut)
                
                linker_alignments[seq_id] = {"qseq": qseq[cut_site:], "hseq": hseq[cut_site:], "midline": midline[cut_site:]}
                LOV2_alignments[seq_id] = {"qseq": qseq[read_out_of_frame:cut_site], "hseq": hseq[read_out_of_frame:cut_site], "midline": midline[read_out_of_frame:cut_site]}

            else:    
                linker_alignments[seq_id] = {"qseq": qseq[:cut_site], "hseq": hseq[:cut_site], "midline": midline[:cut_site]}
                LOV2_alignments[seq_id] = {"qseq": qseq[cut_site:], "hseq": hseq[cut_site:], "midline": midline[cut_site:]}
        else:
            LOV2_start_indel_count +=1

    print(LOV2_start_indel_count, "sequences are excluded, since LOV2 start site could not be found in the ref (due to '-' i.e. insertions at the start of LOV2)")

    return(linker_alignments, LOV2_alignments)
        


def characterize_DMS_blast_alignment(DMS_alignments, ref, data_type = "AA", read_dir = "R1", exclude_not_covered_regions = True):
    """
    Function to characterize the DMS alignments, by counting the number of insertions, deletions and substitutions per position

    args:
    DMS_alignments: dict, with the sequences of insert that is mutated 
    ref: str, reference DNA sequence 
    data_type: str, "AA" or "DNA"
    read_dir: str, "R1" or "R2"

    returns:
    all_variants: dict, with the counts of the variants per position
    indels: pd.DataFrame, with the counts of insertions and deletions per position
    enrichment_counts: pd.DataFrame, with the counts of the variants per position, with the reference sequence masked
    enrichment_relative: pd.DataFrame, with the relative counts of the variants per position, with the reference sequence masked

    """

    all_variants = {}
    seq_with_off_target_indels = 0
    included_seq = 0
    indels = pd.DataFrame(columns = range(len(ref)), index = ["insertion", "deletion"], data = 0)

    ref_prot = translate_dna2aa(ref)

    if data_type == "AA":
        for idx in range(len(ref_prot)):
                all_variants[idx] = {'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 
                                    'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 
                                    'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 
                                    'W':0, 'Y':0, '*':0} ## X: are triplets with missing nucleotides (i.e. in the aligned seq, "-" is present), that would lead to a frameshift
    else: 
        for idx in range(len(ref)):
                all_variants[idx] = {'A':0, 'C':0, 'G':0, 'T':0}
            
    for x in DMS_alignments.values():
        qseq = x["qseq"]
        hseq = x["hseq"]

        if "-" in hseq or "-" in qseq:
            seq_with_off_target_indels += 1
            shift = 0 ## here, we count the shift of the position compared to the reference, that occurs if there is an insertion in the qseq 

            for idx,nt in enumerate(qseq): 
                pos = idx - shift # adjust for the shift in the index, due to prior insertions

                if read_dir == "R2": 
                    pos = len(ref) - (len(hseq)-qseq.count("-")) + pos

                if hseq[idx] == "-":
                    indels.loc["deletion", pos] += 1

                if nt == "-":
                    indels.loc["insertion", pos] += 1
                    shift += 1 ## to correct for the shift in the index, due to the insertion
                
            continue
        #filtered_hseq = "".join([s2 for s1, s2 in zip(qseq, hseq) if s1 != "-"]) ## filter out gaps that are present in the hseq, to exclude insertions that lead to frameshifts (seq errors or off target insertions)

        if data_type == "AA":
            hseq = translate_dna2aa(hseq)

        #print(len(hseq)/3)
        if read_dir == "R2": 
            hseq = hseq[::-1]
            for idx, variant in enumerate(hseq):
                all_variants[len(all_variants)-idx-1][variant] += 1
            included_seq +=1
        else: 
            for idx, variant in enumerate(hseq): 
                all_variants[idx][variant] += 1
            included_seq +=1
       
    indels_freq = indels/len(DMS_alignments)
    print(seq_with_off_target_indels, "sequences with off target indels are excluded")
    print(included_seq, "sequences are included in the enrichment analysis")

    enrichment_df = pd.DataFrame.from_dict(all_variants)

    if exclude_not_covered_regions: 
        enrichment_df = enrichment_df.loc[:,enrichment_df.sum() > 0]

    enrichment_counts, enrichment_relative = mask_ref_in_variants_df(variant_df=enrichment_df, ref_seq=ref if data_type=="DNA" else ref_prot, data_type=data_type, reverse = True if read_dir == "R2" else False)
    

    return all_variants, indels_freq, enrichment_counts, enrichment_relative



def get_linker_variants_from_blast_alignment(linker_alignments, wt_linker = "SG", read_dir = "R1"):
    """
    Function to characterize linker variants from the blast alignments

    args: 
    linker_alignments: dict, with the blast algined sequences (hseq, qseq) of the linker
    wt_linker: str, the wildtype linker AA sequence 
    read_dir: str, "R1" or "R2"

    returns:


    """

    frameshifts = 0
    linkers = {}

    for x in linker_alignments.values():
        qseq = x["qseq"]
        hseq = x["hseq"]

        ##### Exclude frameshift reads 
        if qseq.count("-") % 3 != 0 or hseq.count("-") % 3 != 0:  
            # Insertions (shown as "-" in ref) not multiple of three lead to frameshifts -> exclude these reads
            frameshifts += 1

        ##### WT sequences
        elif qseq == hseq:  # WT linkers with differences in the rest of the sequence are taken into account below
            linkers["wt"] = linkers.get("wt", 0) + 1

        ##### Reads with deletions 
        elif hseq.count("-") > 0 and hseq.count("-") % 3 == 0:  # Deletions that are multiple of 3, not leading to frameshifts
            if hseq.count("-") == 3:  
                # 3 deletions represent substitution of SG linker by a single AA, here we might also get some noise, due to sequencing errors that lead to deletions of 3 Nts but this is acceptable
                hseq_filt = re.sub("-", "", hseq)
                if read_dir == "R2": 
                    linker = translate_dna2aa(hseq_filt)[:len(wt_linker)-1]  # Linker shortened by 3 Nts = 1 AA
                else: 
                    linker = translate_dna2aa(hseq_filt)[-len(wt_linker)+1:]  # Linker shortened by 3 Nts = 1 AA

                linkers[linker] = linkers.get(linker, 0) + 1

            else:  # Deletions in the linker region
                del_count = hseq.count("-") 
                delname = "del-" + str(del_count)
                linkers[delname] = linkers.get(delname, 0) + 1

        ###### Reads with substitutions 
        elif qseq.count("-") == 0:  # Linker was substituted, but no deletions or insertions present
            if read_dir == "R2":
                if hseq[:len(wt_linker) * 3] == qseq[:len(wt_linker) * 3]:
                    linkers["wt"] = linkers.get("wt", 0) + 1
                else: 
                    linker = translate_dna2aa(hseq)[:len(wt_linker)]
                    linkers[linker] = linkers.get(linker, 0) + 1
            else: 
                if hseq[-len(wt_linker) * 3:] == qseq[-len(wt_linker) * 3:]:  
                    # WT linker (but differences in the rest (e.g beginning) of the sequence, thus these did not meet the first criterion)
                    linkers["wt"] = linkers.get("wt", 0) + 1
                else: 
                    linker = translate_dna2aa(hseq)[-len(wt_linker):]
                    linkers[linker] = linkers.get(linker, 0) + 1

        ###### Reads with insertions
        elif qseq.count("-") > 0:
            insertion_len = qseq.count("-") // 3  # AA level
            if read_dir == "R2": 
                linker = translate_dna2aa(hseq)[:len(wt_linker) + insertion_len]
                linkers[linker] = linkers.get(linker, 0) + 1
            else: 
                linker = translate_dna2aa(hseq)[-len(wt_linker) - insertion_len:]
                linkers[linker] = linkers.get(linker, 0) + 1

        #### All other reads
        else:
            print("sequence", hseq, "does not meet any criteria")

    print(frameshifts, "reads excluded due to frameshifts")

    return linkers