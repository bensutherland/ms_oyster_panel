#!/bin/bash
# Uses various .csv files to create a bed file

# Global variables
OUTPUT_FOLDER="04_extract_loci"
INPUT="vcf_selection.csv"
OUTPUT="vcf_selection.bed"

# Prepare bed file
echo "Subtracting 201 bp from the front and adding 200 bp to the back"
echo "...since preparing for a 0-state count"

awk -F, '{ print $1, $2-201, "\t", $2+200 }' $OUTPUT_FOLDER/$INPUT > $OUTPUT_FOLDER/$OUTPUT

