#!/bin/bash

barcode="barcode09"
# Define input variables
raw_fastq_input_folder="/var/lib/minknow/data/basecalling/pass/$barcode"
output_folder="/home/student/anna/DMS_analysis/data/Nanopore/$barcode/filtered_Q15_maxminlen"
output_plot_folder="/home/student/anna/DMS_analysis/output/Nanopore/$barcode/filtered_Q15_maxminlen/quality_control"
reference_file="/home/student/anna/DMS_analysis/data/Nanopore/AraC_S170_LOV_R5_ref.fa"

# Run Python scripts sequentially
echo "################## Running filtering... ##################"
python3 Nanopore_read_filtering.py "$raw_fastq_input_folder" "$output_folder/filtered_fastq"

echo "################## Running alignment... ##################"
python3 Nanopore_alignment.py "$output_folder/filtered_fastq" "$output_folder/minimap2_alignment" "$reference_file"


echo "################## Plotting quality plots... ##################"
python3 Nanopore_quality_control.py "$output_folder/minimap2_alignment" "$output_plot_folder"

echo "################## Running read processing... ##################"
python3 process_Nanopore_reads.py "$output_folder/minimap2_alignment" "$output_folder/processed_reads" "$reference_file"


echo "Pipeline finished!"
