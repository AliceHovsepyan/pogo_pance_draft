
# Illumina
## Final pipeline
### input 
as input, we provide: 
1. a **json** file, in which arguments for filtering, and demultiplexing are specified (thereby, we can later also track which settings were used)

2. fastq files (R1, R2 reads)
#### preprare data:
all of the steps below, are included in the `filter_and_demultiplex_reads.py` file (steps 1-5) and the `analyze_mutation_enrichment.py` or `characterize_indels_from_blast.py` file (step 6)
1. **read filtering:** 	
	- filtering reads: cut reads at first position with error rate > 1%
	- **no** filtering for n_mut, no filtering for read_len
	- we **filter for idx swapping** by just keeping reads with same BC in R1, R2
	- we cut the BC sequence and primer start, so that the reads start (R1) or end (R2) in frame
	- demultiplexing 
	
4. **save demultiplexed reads as fasta files:**
	- for each BC and each Section, reads are saved as **fasta** (in "preprocessed" folder)
	- !!!! importantly, reverse (R2) reads are saved as reverse complented sequences !!!!, thereby, they can be mapped directly on the Plus strand of the same reference sequence, and this simplifies further analysis
	- the reference sequence, per section and barcode is also saved as fasta (in "reference" folder)
	
5. **run blast:**
		
	- we do this on a nucleotide level, since then frameshift reads can still be aligned:
	- below, the bast arguments are shown, but this is all automated within the python script
			
		- **create database with:**
			
			`makeblastdb -in ../../fastq/R36/R36_BC3_Nt_filt_R2_001.fasta -dbtype nucl -out R36_BC3_Nt_filt_R2_001`
			
		- **run blast with:**
			
			`blastn -db R35_BC1_Nt_filt_R2_001 -query ../../fastq/R35/R35_Nt_ref.fasta -out ../../blastoutput/R35_BC1_Nt_filt_R2_001.out -outfmt 15 -max_target_seqs 40000`
			
6. **Analysis of blast reads:**
		
	- filter for reads that are aligned on the LOV2 site i.e. AraC variants, that dont have LOV2 inserted, are here excluded
	- depending on the analysis: 
		**for RL1 (also if combined with RL8)**  run: `characterize_indels_from_blast.py` 
		- separate reads at the LOV2 site in
			1. reads that include the linker sequence (read start until LOV2 (LATTLER))
			2. reads that include LOV2 (read from the start of LOV2 until end (linker not included here))
			
		-> do the analysis of LOV2 and linker reads separately: 
		- **for LOV2 reads**:
			- filter reads with indels (i.e. reads contain “-” in qseq (inseriton) or hseq(deletion)) → frameshift reads
			- store indels per position 
			- calculate mutation enrichment per positon 
		- **for linker reads**
			- **filtering:** reads that have deletions or insertions not multiple of 3 are excluded (frameshifts)
			- then, the linker variants are determined 
			- !!!! what if 3 deletions and 1 insertion????????????????????????? -> fix this! 

		**for RL8/DP6** run  `analyze_mutation_enrichment.py`
		- reads are not separated, since we only induce mutations 
		- exclude reads that are shifted, i.e. do not start at the amplicon site (50-100 reads approx)
		- exclude reads that have insertions or deletions, thereby saving where these indels occur
		- calculate enrichment for the reads per position


## pymol

start pymol in terminal: 
`pymol`


# Nanopore
## Preprocessing and quality assessment

1. Nanopore sequencing with Minknow software
	- thereby, reads are quality filtered so that only reads with a phred score >= 8 are kept
	
2. perform .bam (alignment) files with the Minknow/Nanopore  software
	-> reads were aligned to the respective AraC-S170-LOV2 reference (R2/R5)
	-> bc01-05: R2
	-> bc06-09: R5
	
3. put reads in right frame, i.e.: 
	- if there is a deletion in a read, this is shown as "-"
	- if there is an insertion in a read, this base is skipped
	- the remaining read is kept as it is 

	- also: reads are cut so that all start at the same position (here, we consider reads that include ref pos 9 (arbitrary choice))

	- therefore, we run the `process_Nanopore_reads.py ` file e.g.:
		`python process_Nanopore_reads.py /var/lib/minknow/data/basecalling/pass/barcode05/alignment/ /home/student/anna/DMS_analysis/data/Nanopore/barcode05 /var/lib/minknow/data/AraC_S170_LOV_R2_ref.fa`
	
4. run general quality control plotting using `NanoPlot`
	- either, to run on one (or more) bam files (thereby specifying the files themselves)
		`NanoPlot -t 2 --bam alignment1.bam alignment2.bam alignment3.bam -o bamplots_downsampled`
	- or, to run the analysis on **all** .bam files within a folder using the `Nanopore_quality_control.py` file: 
		`python Nanopore_quality_control.py /var/lib/minknow/data/basecalling/pass/barcode09/alignment /home/student/anna/DMS_analysis/output/Nanopore/barcode09/quality_control`

