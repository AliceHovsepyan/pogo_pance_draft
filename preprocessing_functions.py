
import os
from Bio.SeqIO import QualityIO
import numpy as np
from utils import dna_rev_comp, translate_dna2aa
import pandas as pd


def read_sequences(variant, 
                   catch_left, 
                   catch_right, 
                   base_dir = None, 
                   cutoff_a_read = None, 
                   cutoff_b_read = None, 
                   quality_score = ['!', '"', '#', '$', '%', '&', "'", '(', ')', '*','+', ',', '-', '.', '/', '0', '1', '2', '3', '4', '5'], 
                   return_qualities_ids = False):
    """
    read sequences from fastq files while filtering for quality score (read is aborted at first base with higher error rate than (defaut) 1%)

    variant: variant name of the fastq files, which follow this structure {variant}_R1_001.fastq and {variant}_R2_001.fastq
    catch_left, catch_right: start (end) of the sequence in the forward read (R1) (reverse read (R2)), e.g. Barcodes and/or primers (will not be included in the analysis)
    base_dir: directory where the fastq files are stored (default: current working directory)
    cutoff_a_read: at which position to cut off all forward reads that already went through quality score filtering (= max length of R1 reads) (default: None)
    cutoff_b_read: at which position to cut off all backward reads that already went through quality score filtering (= max length of R2 reads) (default: None)
    return_qualities_ids =  whether or not to return lists of R1 qualities, R2 qualities, R1 ids, R2 ids (default: False)
    quality_score: list of quality scores, at which the reads should be aborted (default: 1% error rate)

    returns: list of  R1 sequences, R2 sequences (and optionally R1 qualities, R2 qualities, R1 ids, R2 ids)
    """

    if not base_dir:
        base_dir = os.getcwd() + "/data/fastq/"

    a_sequences = []
    b_sequences = []
    a_qualities = []
    b_qualities = []
    a_ids = []
    b_ids = []

    with open(f'{base_dir}/{variant}_R1_001.fastq', "rt") as a_file, open(f'{base_dir}/{variant}_R2_001.fastq', "rt") as b_file:

        a_reader = QualityIO.FastqGeneralIterator(a_file)
        b_reader = QualityIO.FastqGeneralIterator(b_file)
        
        for total_read, (a, b) in enumerate(zip(a_reader, b_reader)):
                
                a_id, a_seq, a_qual = a
                b_id, b_seq, b_qual = b
                cutoff_a = find_(a_qual, quality_score)
                cutoff_b = find_(b_qual, quality_score)

                if cutoff_a_read and catch_left in a_seq: # cut off a_seq to an (arbitrary) chosen maximum length (=cutoff_a_read) after the catch_left sequence
                    if cutoff_a > (a_seq.index(catch_left) + cutoff_a_read):
                        cutoff_a = a_seq.index(catch_left)  + len(catch_left) + cutoff_a_read 
                
                if cutoff_b_read and dna_rev_comp(catch_right) in b_seq: 
                    if cutoff_b > (b_seq.index(dna_rev_comp(catch_right)) + cutoff_b_read):
                        cutoff_b = b_seq.index(dna_rev_comp(catch_right))+ len(catch_right) + cutoff_b_read

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


def demultiplex_reads(a_seqs:list, 
                      b_seqs:list,
                      Barcodes:dict, 
                      Primer_seq:dict, 
                      Primer_out_of_frame:dict,
                      used_Barcodes:list, 
                      Sections:list, 
                      max_mismatch_primerseq:int = 5, 
                      a_ids:list = None, 
                      b_ids:list = None, 
                      cut_BC_seq = True,
                      cut_primer_start=True,
                      catch_left = "",
                      catch_right = "",
                      include_only_complete_reads = False):
    """
    demultiplex reads from fastq-files, if different samples were pooled and the region of interest was divided into sections for sequencing
    
    a_seqs, b_seqs: list of forward reads (R1) and reverse reads (R2)
    ref_gene: reference DNA sequence
    Barcodes: dictionary with the forward and reverse Barcode sequences, following the structure {BC1_fwd : seq, BC1_rev : seq, BC2_fwd : seq, ... }
    Primer_seq: dictionary with the fwd and rev primer sequences for each section, following the structure {S1_fwd : seq, S1_rev : seq, S2_fwd : seq, ... }
    Primer_out_of_frame: dictionary with the number of nucleotides at the beginning of the primer seq before a triplet starts, following the structure {S1_fwd : int, S1_rev : int, S2_fwd : int, ... } 
    used_Barcodes: list of Barcodes from the Barcodes dictionary that were used for the sequencing (should match with names in the Barcodes dict, e.g. BC1, BC2, ...)
    Sections: list of sections that were sequenced (should match with names in the Primer_seq dict, e.g. S1, S2, ...)
    max_mismatch_primerseq: maximum number of mismatches allowed in the primer sequences (default: 5) (to keep reads that contain mutations in the primer seq)
    a_ids, b_ids: list of ids for the forward and reverse reads (default: None), if None, no ids are returned
    cut_BC_seq: whether or not to cut the BC seq from the reads 
    cut_primer_start: used if cut_BC_seq=True, then, if True, the Nt number specified in Primer_out_of_triplets is additionally cut from the sequence start, to keep the reads in frame
    catch_left, catch_right: start (end) of the sequence in the forward read (R1) (reverse read (R2)), after which the sequence should be cut (will not be included in the analysis)
    (default: "", i.e. no cutting)

    returns: dictionary with the reads for each sample and section, optionally also the ids
    """

    read_Dict = {}
    ids_Dict = {}

    ## split the reads into the samples according to Barcode and Section -> thereby keeping the forward and reverse reads together
    for Barcode in used_Barcodes: 

        for Section in Sections:

            fwd_BC_Primer_seq = Barcodes[Barcode + "_fwd"] + Primer_seq[Section+"_fwd"]
            rev_BC_Primer_seq = Barcodes[Barcode + "_rev"] + Primer_seq[Section+"_rev"] 

            ### select the reads that contain the forward and reverse BC + primer sequences, thereby allowing for n mismatches in the primer sequences but no errors in BCs
            fwd_idxs = []
            rev_idxs = []

            for a_idx, seq in enumerate(a_seqs):
                a_mismatch_to_primer_seq = sum([sequence!=primer_ref for sequence, primer_ref in zip(seq[len(Barcodes[Barcode + "_fwd"]):len(fwd_BC_Primer_seq)], Primer_seq[Section+"_fwd"])])
                if seq[:len(Barcodes[Barcode + "_fwd"])] == Barcodes[Barcode + "_fwd"] and a_mismatch_to_primer_seq <= max_mismatch_primerseq:
                    fwd_idxs.append(a_idx)

            for b_idx, seq in enumerate(b_seqs):
                b_mismatch_to_primer_seq = sum([sequence!=primer_ref for sequence, primer_ref in zip(seq[len(Barcodes[Barcode + "_rev"]):len(rev_BC_Primer_seq)], Primer_seq[Section+"_rev"])])
                if seq[:len(Barcodes[Barcode + "_rev"])] == Barcodes[Barcode + "_rev"] and b_mismatch_to_primer_seq <= max_mismatch_primerseq:
                    rev_idxs.append(b_idx)
            
            indexes = set(
                [idx for idx in fwd_idxs if b_seqs[idx][:len(Barcodes[Barcode + "_rev"])] == Barcodes[Barcode + "_rev"]]  +  
                [idx for idx in rev_idxs if a_seqs[idx][:len(Barcodes[Barcode + "_fwd"])] == Barcodes[Barcode + "_fwd"]])## only keep reads that match in the fwd and rev BC seqs
                
            print(sum([len(b_seqs[fwd_i]) < len(Barcodes[Barcode + "_rev"]) for fwd_i in fwd_idxs]), "b reads are empty") ## reads that are only in the reverse list
            print(sum([len(a_seqs[rev_i]) < len(Barcodes[Barcode + "_fwd"]) for rev_i in rev_idxs]), "a reads are empty") ## reads that are only in the reverse list

            print(len(indexes), "reads with matching BC and primer seq")
            print(len(set(fwd_idxs+ rev_idxs)) - len(indexes), "reads with index swapping")

            a_seq_Bc_Sec = [a_seqs[i] for i in indexes]
            b_seq_Bc_Sec = [b_seqs[i] for i in indexes]

            print(Barcode, Section, len(a_seq_Bc_Sec), "total fwd reads")

            if a_ids and b_ids:
                a_ids_Bc_Sec = [a_ids[i].split(" ")[0] for i in indexes]
                b_ids_Bc_Sec = [b_ids[i].split(" ")[0]  for i in indexes]

    
            if cut_BC_seq: 
                cutoff_a = len(Barcodes[Barcode + "_fwd"]) if not cut_primer_start else len(Barcodes[Barcode + "_fwd"]) + Primer_out_of_frame[Section + "_fwd"]
                cutoff_b = len(Barcodes[Barcode + "_rev"]) if not cut_primer_start else len(Barcodes[Barcode + "_rev"]) + Primer_out_of_frame[Section + "_rev"]

                a_seq_Bc_Sec = [a[cutoff_a:] if len(a)>=cutoff_a else "" for a in a_seq_Bc_Sec]
                b_seq_Bc_Sec = [b[cutoff_b:] if len(b)>=cutoff_b else "" for b in b_seq_Bc_Sec]
            
            ## cut sequences at the catch_left and catch_right positions
            if include_only_complete_reads: ## only include reads that contain the full sequence (i.e. catch_left **and** catch_right is present)
                a_seq_Bc_Sec = [read[read.index(catch_left)+len(catch_left):read.index(catch_right)] if catch_left in read and catch_right in read else "" for read in a_seq_Bc_Sec ]

                b_seq_Bc_Sec = [read[read.index(dna_rev_comp(catch_right))+len(catch_right):read.index(dna_rev_comp(catch_left))] if dna_rev_comp(catch_right) in read  and dna_rev_comp(catch_left) in read else "" for read in b_seq_Bc_Sec ]

            else: ## cut sequences at the catch_left and catch_right positions, reads do not have to be complete 
                a_seq_Bc_Sec = [a[a.index(catch_left)+len(catch_left):] if catch_left in a else "" for a in a_seq_Bc_Sec]
                b_seq_Bc_Sec = [b[b.index(dna_rev_comp(catch_right))+len(catch_right):] if dna_rev_comp(catch_right) in b else "" for b in b_seq_Bc_Sec]

            read_Dict[f"{Barcode}_{Section}_R1"] = a_seq_Bc_Sec
            read_Dict[f"{Barcode}_{Section}_R2"] = b_seq_Bc_Sec

            if a_ids and b_ids:
                ids_Dict[f"{Barcode}_{Section}_R1"] = a_ids_Bc_Sec
                ids_Dict[f"{Barcode}_{Section}_R2"] = b_ids_Bc_Sec

            print(f"################# Completed {Barcode} {Section} #################")

        print(f"################# Completed {Barcode} #################")
    if a_ids and b_ids:
        return  read_Dict, ids_Dict
    else:
        return read_Dict 



def find_reference_seq(ref_gene, 
                       Primer_seq, 
                       Section, 
                       Primer_out_of_frame):
    """
    find the reference sequence for a given section within a reference gene, based on the primer sequences (takes into account that the primers can be out of frame)

    ref_gene = reference gene sequence
    Primer_seq: dictionary with the fwd and rev primer sequences for each section, following the structure {S1_fwd : seq, S1_rev : seq, S2_fwd : seq, ... }
    Section: Section of interest
    Primer_out_of_frame = dictionary with the number of nucleotides at the beginning of the primer seq before a triplet starts, following the structure {S1_fwd : int, S1_rev : int, S2_fwd : int, ... }

    returns: reference sequence for the given section
    """ 
    tripl_st = Primer_out_of_frame[Section+"_fwd"]
    tripl_end = Primer_out_of_frame[Section+"_rev"]
    primer_fwd = Primer_seq[Section + "_fwd"][tripl_st:]
    primer_rev = dna_rev_comp(Primer_seq[Section+"_rev"][tripl_end:])
    
    ref_gene_section = ref_gene[ref_gene.index(primer_fwd):ref_gene.index(primer_rev)+len(primer_rev)]

    return ref_gene_section

