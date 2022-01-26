#!/bin/bash
# Use bed file to extract relevant sections of a fasta

# Global variables
INPUT_FOLDER="02_input_data"
REF_GENOME="GCA_000297895.1_oyster_v9_genomic.fna"
OUTPUT_FOLDER="04_extract_loci"
BED_FILE="vcf_selection.bed"
OUTPUT="vcf_selection.fa"

# Remove any existing output
echo "Removing $OUTPUT to make new file" 
rm $OUTPUT_FOLDER/$OUTPUT

# Extract from fasta
bedtools getfasta -fi $INPUT_FOLDER/$REF_GENOME -bed $OUTPUT_FOLDER/$BED_FILE > $OUTPUT_FOLDER/$OUTPUT


