#!/bin/bash

# Define input and output directories
INPUT_DIR="/var/lib/minknow/data/basecalling/pass/barcode01/alignment"
OUTPUT_DIR="/home/student/anna/DMS_analysis/data/Nanopore"

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through all BAM files in the input directory
for file in "$INPUT_DIR"/*.bam; do
    if [[ -f "$file" ]]; then  # Check if file exists to avoid errors
        filename=$(basename "$file" .bam)  # Extract filename without extension
        samtools view -h -o "$OUTPUT_DIR/${filename}.sam" "$file"
        echo "Converted: $file -> $OUTPUT_DIR/${filename}.sam"
    fi
done

echo "All BAM files have been converted to SAM format."