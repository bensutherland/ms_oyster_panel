#!/bin/bash
# Uses various .csv files to create a bed file

# Global variables
OUTPUT_FOLDER="04_extract_loci"
INPUT="vcf_selection.csv"
OUTPUT="vcf_selection.bed"

# Prepare bed file
echo "Subtracting 201 bp from the front and adding 200 bp to the back"
echo "...since preparing for a 0-state count"

# Keep the chromosome (1), the lower range, upper range, and marker info
awk -F, '{ print $1 "\t" $2-201 "\t" $2+200 "\t" "mname_"$3 }' $OUTPUT_FOLDER/$INPUT |
    
    # Remove extra information in last field trailing the colon
    awk -F: '{ print $1 }' - >  $OUTPUT_FOLDER/$OUTPUT

