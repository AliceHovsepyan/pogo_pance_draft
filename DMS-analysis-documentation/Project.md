### RL8
#### Filtering 
- reads that differ at more than 10 Nt positions are excluded (`filter_for_n_mut` )
	--> these are likely frameshift reads, i.e. sequencing errors (bases are skipped, or read twice) leading to a high number of "mutations" 
- reads with len < 15 are excluded (`filter_for_read_len` )
	--> these pass the `filter_for_n_mut` filter, since these are short enough, resulting in high "mutation rates" at the sequence start
- we abort reads at the first base with a error rate > 1% (`quality_score`)
- we filter reads that lack the *As*LOV2 domain, by looking at their sequences 

**for multiplexed sequencing experiments**
- we filter for reads that contain the correct barcode in the fwd and rev read (we also exclude reads for which the fwd or rev read is empty)
	--> to filter out reads with index swapping

#### Analysis
##### of positional effects
- there are regions with higher/lower mut rate --> why??
	- hybridization energy of our retron sequence (only the part that hybridizes) with the target sequence? 
		- calculate delta G with UnaFold (online tool)
		- see:  [UnaFold Server](https://www.unafold.org/Dinamelt/applications/two-state-melting-hybridization.php#)
		![](https://lh7-rt.googleusercontent.com/slidesz/AGV_vUeAxhvzL6V6CKu9kUH2PP32Asj5JRqqLisJw_dNxQpj5T3Vynxjun3b3xcjBvTJiN61o5HDdGZ5uFFENwrLJL8njZWIDTk0QAeRLxFGU1sbwDX7YAk8lL3tJPjVEx7NZPKQSzsCSw=s2048?key=EPQFgKB0w4R2jR544RnKdg)
		→ less negative ∆G = less stable binding 
		→ only the retron region that hybridizes to the target seq is considered
		→ does not really correlate… 
	
	- Are low mutation rates correlated with secondary DNA structures within the retrons? (ViennaRNAfold)
			-->  hypothesis:  high ∆G → strong secondary structures within retrons → impairs binding to the target sequence → lower mutation rates 
			-->  here, we consider only the retron region that binds to the target sequence
			-->  but we can even see an anti correlation → hypothesis is not confirmed
			-->  package documentation: [https://www.tbi.univie.ac.at/RNA/RNAfold.1.html#heading11](https://www.tbi.univie.ac.at/RNA/RNAfold.1.html#heading11) 
			function call: 
		`RNAfold -i AraC_S170_LOV_DMS_Retron_Pool.fasta --noconv --param=DNA -oviennaRNAfoldAraC_S170_LOV_DMS_Retron_Pool


### DP6
- same filtering as RL8

### RL1
- we abort reads at the first base with a error rate > 1% (`quality_score`)
- we do not `filter_for_n_mut` and `filter_for_read_len`
- instead, we filter for reads containing our region of interest 
##### CRISPResso
- we can also use crispresso to look for indels 
- thereby, we consider the linker (indel target site) as "CRISPR cut site"
- **preparation of rev read** (reverse complement the sequences in the fastq files to get the correct orientation of the sequences (in command line))
	`seqkit seq -r -p test.fastq -o revcomp_test.fastq -t DNA
**Function calls:** 
-   for the forward read only (only first 150 positions of the LOV gene (+beginning of amplicon))
	`CRISPResso --fastq_r1 RL1_R1_001.fastq --amplicon_seq CGCCGCATGGAAGCGATTAACGAAAGCAGCGGTTTAGCCACAACGCTGGAACGCATTGAAAAGAATTTCGTAATCACAGACCCGCGCCTTCCCGACAATCCAATTATTTTTGCGTCCGATAGCTTCCTGCAATTAACCGAATACAGCCGCGAAGAAATTCTGGGTCGTAATTGTCGCTTCCTT -n RL1-R1_crispresso_result_final --guide_seq GAAGCGATTAACGAAAGCAGCGGT --cleavage_offset 0 --amplicon_min_alignment_score 50 --min_average_read_quality 20 --min_bp_quality_or_N 10 --offset_around_cut_to_plot 25 --window_around_sgrna 10 --min_frequency_alleles_around_cut_to_plot 0.00005 --needleman_wunsch_gap_incentive -20 --needleman_wunsch_gap_extend 0

- for both (R1 and R2), however due to low coverage, many reads are not aligned well enough and thus excluded from the analysis
	`CRISPResso --fastq_r1 RL1_R1_001.fastq --fastq_r2 RL1_R2_001.fastq --amplicon_seq CGCCGCATGGAAGCGATTAACGAAAGCAGCGGTTTAGCCACAACGCTGGAACGCATTGAAAAGAATTTCGTAATCACAGACCCGCGCCTTCCCGACAATCCAATTATTTTTGCGTCCGATAGCTTCCTGCAATTAACCGAATACAGCCGCGAAGAAATTCTGGGTCGTAATTGTCGCTTCCTTCAGGGGCCAGAGACTGACCGTGCTACGGTACGCAAAATCCGCGACGCAATCGACAATCAAACGGAAGTCACGGTTCAGTTGATTAACTATACGAAGAGCGGAAAAAAATTCTGGAATTTATTTCACTTGCAGCCTATGCGTGACCAGAAGGGCGATGTCCAGTATTTCATTGGCGTTCAGCTTGATGGTACCGAGCATGTTCGCGATGCTGCGGAGCGTGAAGGTGTAATGTTAATTAAAAAGACTGCTGAAAACATTGATGAAGCGGCCAAAGGGAGCCTGCATCCGCCGATGGATAACCGCGTG -n RL1_crispresso_result --guide_seq GAAGCGATTAACGAAAGCAGCGGT,TGAAAACATTGATGAAGCGGCCAAA --cleavage_off





### POGO 

#### Filtering 
- reads that differ at more than 10 Nt positions are excluded (`filter_for_n_mut` )
	--> these are likely frameshift reads, i.e. sequencing errors (bases are skipped, or read twice) leading to a high number of "mutations" 
- Primer sequences are cut from the sequences (thus, we do not need `filter_for_read_len`) 
- we abort reads at the first base with a error rate > 1% (`quality_score`)
- we filter for reads that contain the correct barcode in the fwd and rev read (we also exclude reads for which the fwd or rev read is empty)
	--> to filter out reads with index swapping

