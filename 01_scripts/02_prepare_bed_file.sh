#!/bin/bash
# Uses various .csv files to create a bed file

# Global variables
RAW_FOLDER="02_input_data"
INPUT_VCF="populations.snps.vcf"
MARKER_FOLDER="04_extract_loci"
OUTPUT="vcf_selection.csv"

# Clearing space
echo "Removing $OUTPUT to make new file" 
rm $MARKER_FOLDER/$OUTPUT

# Prepare bed file
cat $MARKER_FOLDER/top*.csv | 
    
    # Troublshooting
    #head -n 5 |
 
    # Keep only the first column for the mname and remove duplicates
    awk -F, '{ print $1 }' - | 
    sort |
    uniq |

    # Remove additional details from the mname
    awk -F_ '{ print $1":" }' - | 

    # For each mname, extract the relevant section of the vcf
    while read i
    do
        echo $i

        grep $i $RAW_FOLDER/$INPUT_VCF |
        awk '{ print $1 "," $2 "," $3 "," $4 "," $5 }' - >> $MARKER_FOLDER/$OUTPUT
 
    done    

