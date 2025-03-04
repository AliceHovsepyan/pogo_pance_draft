#!/bin/bash

# Define arguments
filepath="/home/student/anna/DMS_analysis/data/fastq/P0111_RetronLib" ## path to the fastq files

# Run the Python script with arguments
echo "Running Python script filter_and_demultiplex_reads.py:"
python3 filter_and_demultiplex_reads.py "$filepath" 