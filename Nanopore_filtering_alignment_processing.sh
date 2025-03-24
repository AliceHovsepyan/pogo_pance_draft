#!/bin/bash

for barcodenumber in 06 07 08 09 
do 
    barcode="barcode"$barcodenumber
    # Define input variables
    raw_fastq_input_folder="/var/lib/minknow/data/basecalling/pass/$barcode" #"/var/lib/minknow/data/P0115/basecalling/pass/$barcode"
    output_folder="/home/student/anna/DMS_analysis/data/Nanopore_P0109/$barcode/highly_accurate_basecalling/filtered_Q20_maxminlen"
    output_plot_folder="/home/student/anna/DMS_analysis/final_output/Nanopore_P0109/$barcode/highly_accurate_basecalling/filtered_Q20_maxminlen/quality_control"
    reference_file="/home/student/anna/DMS_analysis/data/Nanopore_P0109/AraC_S170_LOV_R5_ref.fa"

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
done