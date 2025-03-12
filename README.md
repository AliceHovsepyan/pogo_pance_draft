
# Illumina

## Final pipeline

### Input 
as input, we provide: 
1. a **json** file, in which arguments for filtering, and demultiplexing should be specified (thereby, we can later also track which settings were used)

2. fastq files (forward (R1), reverse (R2) reads) -> the filenames should follow the notation: ({variant name}_{read direction}_001.fastq, e.g. DP6_R1_001.fastq)

#### Analysis
All of the steps below, are automatized in the `preprocess_and_align_illumina_reads.py` file (steps 1-3) and the `analyze_mutation_enrichment.py` or `characterize_indels_from_blast.py` file (step 4). 
If you provide the correct input, you can just run these files as follows: 

`python preprocess_and_align_illumina_reads.py input_folder --save_ref`

Please make sure to set the correct parameters at the beginning of the `analyze_mutation_enrichment.py` or `characterize_indels_from_blast.py` files!

Below, you can find some additional information on the steps included in these files:

1. **read filtering:** 	
	- **quality filtering**: cut reads at first position with error rate > **1%**
	- **filtering for idx swapping**: only keep reads with same Barcode (BC) in R1, R2
	- **demultiplexing**: reads are demultiplexed by provided Barcode and Primer sequences
	- we **cut** the Barcode sequence and primer start, so that the reads are in frame
	
2. **save demultiplexed reads as fasta files:**
	- for each Barcode and each Section, reads are saved as **fasta** files (in the "/preprocessed" folder)
	- importantly, **reverse** (R2) reads are saved as **reverse complented** sequences !!!!, thereby, they can be mapped directly on the Plus strand of the reference sequence, which simplifies further analysis
	- if the `--save_ref` tag is provided (**recommended**), the **reference sequence** per section and barcode is also saved as fasta (in "/reference" folder) (based on `amplicon` provided in `config.json` file). If you want to provide the reference sequence yourself, please don't set `--save_ref`, but make sure to provide the reference fasta file(s) within the /reference folder with the right filename(s): {variant}_{Barcode}_{Section}_Nt_filt_ref.fasta (one reference per Barcode and Section (same reference for R1 and R2 reads))
	
3. **run blast:**
		
	- we align the reads on the **nucleotide** level, since then, frameshift reads can still be aligned
	- this includes (below, the bast arguments are shown, but this is all automated within the python script `preprocess_and_align_illumina_reads.py`): 
		1. **creating a blast database** with: 
			
			`makeblastdb -in preprocessed/DP6_BC1_S1_Nt_filt_R1.fasta -dbtype nucl -out blast/db/DP6_BC1_S1_Nt_filt_R1`
			
		2. **run blast with:**
			
			`blastn -db blast/db/DP6_BC1_S1_Nt_filt_R1 -query references/DP6_BC1_S1_Nt_filt_ref.fasta -out blast/alignments/DP6_BC1_S1_Nt_filt_R1.out -outfmt 15 -max_target_seqs 100000`

			thereby, `-outfmt 15` stores the output as json files
			
4. **Analysis of blast reads:**
		
	- please make sure to **set the right parameters** within the file, matching your requirements (e.g. the right filepath, read directory, whether the analysis should be performed on the DNA, AA or Codon level, ...)
	- filter for reads that are aligned on the insert (e.g. LOV2) site i.e. reads that do not span the insertion site, and AraC variants, that dont have an insert, are excluded
	- depending on the analysis: 
		1. **for Retron libary 1 (RL1) (also if combined with RL8)** run: `characterize_indels_from_blast.py` 
			- separate reads at the insert (e.g LOV2) site in
				1. sequences that **include the linker sequence**, i.e. read start until insert after linker (e.g. LOV2: LATTLER)
				2. sequences that include the **insert** (read from the start of the insert until read end (linker **not included** here))
				
			the analysis of the insert and linker sequences is performed separately: 
			- **for sequences spanning the insert**:
				- filter out reads with **indels** (i.e. reads contain “-” in qseq (inseriton) or hseq(deletion)) → frameshift reads, thereby saving where indels occur
				- calculate **mutation enrichment** per positon 
			- **for linker reads**
				- **filtering:** reads that have deletions or insertions not multiple of 3 are excluded (frameshifts) -> we cannot tell whether these frameshifts are due to sequencing errors or biologically accurate -> thus, we just filter out all of them -> these variants would anyways lead to nonsense proteins
				- then, the **linker variants** are determined 

		2. **for RL8/DP6 DMS screens** run `analyze_mutation_enrichment.py` 
			- reads are not separated, since we only induce mutations 
			- exclude reads that are **shifted**, i.e. do not start at the amplicon site (only a small fraction)
			- exclude reads that have **indels**, thereby saving where these indels occur
			- calculate **mutation enrichment** for the reads per position
			- calculate **mutagenic spectrum**, i.e. what mutations are induced


## pymol

start pymol in terminal: 
`pymol`

# Nanopore

### Sequencing
1. Nanopore sequencing with **Minion/Minknow** software
	- thereby, basecalling is performed using the model with the **highest accuracy**

### Read processing
Please run the `Nanopore_filtering_alignment_processing.sh` script, which automatically runs all steps (1-5) below (make sure to set the right paths within the bash script):

2. **Quality filtering** using `chopper` -> filer for read length and average read quality (Q>20)
	
3. **alignment** using `minimap2` 

4. **processing**, i.e. force reads in right frame: 
	- if there is a deletion in a read, this is shown as "-"
	- if there is an insertion in a read, this base is skipped
	- the remaining read is kept as it is 

	- also: reads are cut so that all start at the same position (here, we consider reads that include ref pos 9 (arbitrary choice))

	- to run this processing step seperately, you can run the `process_Nanopore_reads.py ` file e.g.:
		`python process_Nanopore_reads.py /var/lib/minknow/data/basecalling/pass/barcode05/alignment/ /home/student/anna/DMS_analysis/data/Nanopore/barcode05 /var/lib/minknow/data/AraC_S170_LOV_R2_ref.fa`
	
5. perform plotting for **quality control** using `NanoPlot`
	- to run this separately, you can run: 
		- either, to run on one (or more) bam files (thereby specifying the files themselves)
			`NanoPlot -t 2 --bam alignment1.bam alignment2.bam alignment3.bam -o bamplots_downsampled`
		- or, to run the analysis on **all** .bam files within a folder using the `Nanopore_quality_control.py` file: 
			`python Nanopore_quality_control.py /var/lib/minknow/data/basecalling/pass/barcode09/alignment /home/student/anna/DMS_analysis/output/Nanopore/barcode09/quality_control`

